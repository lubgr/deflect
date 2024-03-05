package deflect

import (
	"errors"
	"fmt"
	"log"
	"slices"

	"gonum.org/v1/gonum/floats/scalar"
	"gonum.org/v1/gonum/mat"
)

// NewLinearProblemSolver creates a linear solver for boundary value problems.
func NewLinearProblemSolver() ProblemSolver {
	return &linearSolver{}
}

type linearSolver struct {
	eqn                matrices
	primary, reactions []NodalValue
	dim, constrained   int
	err                error
}

type matrices struct {
	k                  *mat.SymDense
	rhs, d             *mat.VecDense
	k11, k22           mat.Symmetric
	k12                *mat.Dense
	d1, d2, rhs1, rhs2 *mat.VecDense
	scratch            *mat.VecDense
}

func (s *linearSolver) Solve(
	p *Problem,
	layout IndexLayout,
	strategy EquationSolver,
) (primary, reactions []NodalValue, err error) {
	if err := s.initialise(len(layout.indices), len(p.Dirichlet)); err != nil {
		return nil, nil, err
	}

	log.Printf("Start assembly and solution for %v unknowns, %v with BC", s.dim-s.constrained, s.dim)

	k, rhs, d, scratch := s.eqn.k, s.eqn.rhs, s.eqn.d, s.eqn.scratch
	k11, k12, k22 := s.eqn.k11, s.eqn.k12, s.eqn.k22
	d1, d2, rhs1, rhs2 := s.eqn.d1, s.eqn.d2, s.eqn.rhs1, s.eqn.rhs2

	for _, e := range p.Elements {
		s.acc(e.Assemble(layout.indices, k, rhs, d))
	}

	var hasNonZeroDirichlet bool

	for _, bc := range p.Dirichlet {
		hasNonZeroDirichlet = hasNonZeroDirichlet || scalar.EqualWithinAbs(bc.Value, 0, 1e-12)
		s.acc(applyNodalBC(bc, layout.indices, d))
	}

	for _, bc := range p.Neumann {
		s.acc(applyNodalBC(bc, layout.indices, rhs))
	}

	for _, transform := range p.EqTransforms {
		s.acc(transform.Pre(layout.indices, k, rhs, d))
	}

	if s.err != nil {
		// Don't attempt to continue if something went wrong so far
		return nil, nil, s.err
	}

	s.extractK12()

	if hasNonZeroDirichlet {
		// Computes [rhs_2 - k_21 d_1] to account for the known, constrained variables
		scratch.Reset()
		scratch.ReuseAsVec(s.dim - s.constrained)
		scratch.MulVec(k12.T(), d1) // T() doesn't copy, it returns a view
		rhs2.SubVec(rhs2, scratch)
	}

	errSolve := strategy.SolveLinearSystem(k22, rhs2, d2)
	if errSolve != nil {
		return nil, nil, fmt.Errorf("failed to solve assembled linear system: %w", errSolve)
	}

	//  Computes [rhs_1] = [k_11 d_1 + k_12 d_2]:
	scratch.Reset()
	scratch.ReuseAsVec(s.constrained)
	scratch.MulVec(k12, d2)
	rhs1.AddVec(rhs1, scratch)
	scratch.MulVec(k11, d1)
	rhs1.AddVec(rhs1, scratch)

	for _, transform := range p.EqTransforms {
		s.acc(transform.Post(layout.indices, rhs, d))
	}

	// Group primary/secondary solution with symbolic index
	for i := 0; i < s.dim; i++ {
		idx := layout.inverse[i]
		s.primary[i] = NodalValue{Index: idx, Value: d.AtVec(i)}
		s.reactions[i] = NodalValue{Index: idx, Value: rhs.AtVec(i)}
	}

	return s.primary, s.reactions, s.err
}

func (s *linearSolver) initialise(dim, constrained int) error {
	s.err = nil

	if dim == constrained {
		return fmt.Errorf("all %v degrees of freedom have Dirichlet BC, no need to solve this", dim)
	}

	s.dim = dim
	s.constrained = constrained

	s.eqn = formMatrices(dim, constrained, &s.eqn)
	s.primary = slices.Grow(s.primary, dim)[0:dim]
	s.reactions = slices.Grow(s.reactions, dim)[0:dim]

	return nil
}

// formMatrices allocates or re-shapes existing matrices and defines the views over subsets of the
// backing matrices. When possible, existing buffers are re-used.
func formMatrices(dim, constrained int, prior *matrices) matrices {
	eqn := prior

	if eqn.k == nil {
		eqn.k = mat.NewSymDense(dim, nil)
		eqn.rhs = mat.NewVecDense(dim, nil)
		eqn.d = mat.NewVecDense(dim, nil)
		eqn.k12 = mat.NewDense(constrained, dim-constrained, nil)
		eqn.scratch = mat.NewVecDense(max(dim-constrained, constrained), nil)
	}

	for _, matrix := range []mat.Reseter{eqn.k, eqn.rhs, eqn.d, eqn.k12, eqn.scratch} {
		matrix.Reset()
	}

	eqn.k.ReuseAsSym(dim)
	eqn.rhs.ReuseAsVec(dim)
	eqn.d.ReuseAsVec(dim)
	eqn.k12.ReuseAs(constrained, dim-constrained)
	eqn.scratch.ReuseAsVec(max(dim-constrained, constrained))

	// The system k·d = rhs is partitioned as
	//   ⎡k_11 k_12⎤⎡d_1⎤  ⎡rhs_1⎤
	//   ⎣k_21 k_22⎦⎣d_2⎦  ⎣rhs_2⎦
	// where d_1 are all Dirichlet-prescribed values, and d_2 are free primary nodal values to be
	// solved for. We solve these through the linear system
	//   [k_22][d_2] = [rhs_2 - k_21 d_1]
	// and then use the solution vector to compute the reaction forces:
	//   [rhs_1] = [k_11 d_1 + k_12 d_2]
	eqn.k11 = eqn.k.SliceSym(0, constrained)
	eqn.k22 = eqn.k.SliceSym(constrained, dim)

	eqn.d1 = eqn.d.SliceVec(0, constrained).(*mat.VecDense)
	eqn.d2 = eqn.d.SliceVec(constrained, dim).(*mat.VecDense)
	eqn.rhs1 = eqn.rhs.SliceVec(0, constrained).(*mat.VecDense)
	eqn.rhs2 = eqn.rhs.SliceVec(constrained, dim).(*mat.VecDense)

	return *eqn
}

func (s *linearSolver) acc(op error) {
	s.err = errors.Join(s.err, op)
}

func (s *linearSolver) extractK12() {
	for i := 0; i < s.constrained; i++ {
		for j := s.constrained; j < s.dim; j++ {
			s.eqn.k12.Set(i, j-s.constrained, s.eqn.k.At(i, j))
		}
	}
}