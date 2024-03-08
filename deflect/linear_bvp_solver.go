package deflect

import (
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
	eqn              matrices
	d, r             []NodalValue
	dim, constrained int
}

type matrices struct {
	k              *mat.SymDense
	r, d           *mat.VecDense
	k11, k22       mat.Symmetric
	k12            *mat.Dense
	d1, d2, r1, r2 *mat.VecDense
	scratch        *mat.VecDense
}

func (s *linearSolver) Solve(
	p *Problem,
	indices EqLayout,
	strategy EquationSolver,
) ([]NodalValue, []NodalValue, error) {
	if err := s.initialise(indices.EqSize(), len(p.Dirichlet)); err != nil {
		return nil, nil, err
	}

	log.Printf("Start assembly and solution for %v unknowns, %v with BC", s.dim-s.constrained, s.dim)

	k, r, d, scratch := s.eqn.k, s.eqn.r, s.eqn.d, s.eqn.scratch
	k11, k12, k22 := s.eqn.k11, s.eqn.k12, s.eqn.k22
	d1, d2, r1, r2 := s.eqn.d1, s.eqn.d2, s.eqn.r1, s.eqn.r2

	for _, e := range p.Elements {
		e.Assemble(indices, k, r, d)
	}

	var hasNonZeroDirichlet bool

	for _, bc := range p.Dirichlet {
		hasNonZeroDirichlet = hasNonZeroDirichlet || scalar.EqualWithinAbs(bc.Value, 0, 1e-12)
		d.SetVec(indices.MapOne(bc.Index), bc.Value)
	}

	for _, bc := range p.Neumann {
		r.SetVec(indices.MapOne(bc.Index), bc.Value)
	}

	for _, transform := range p.EqTransforms {
		transform.Pre(indices, k, r, d)
	}

	if err := indices.Failure(); err != nil {
		// Don't attempt to continue if something went wrong so far
		return nil, nil, fmt.Errorf("failed to assemble global matrices: %w", err)
	}

	s.extractK12()

	if hasNonZeroDirichlet {
		// This computes [r_2 - k_21 d_1] (see partitioning below) to account for the known, constrained
		// primary nodal values:
		scratch.Reset()
		scratch.ReuseAsVec(s.dim - s.constrained)
		scratch.MulVec(k12.T(), d1) // T() doesn't copy, it returns a view
		r2.SubVec(r2, scratch)
	}

	errSolve := strategy.SolveLinearSystem(k22, r2, d2)
	if errSolve != nil {
		return nil, nil, fmt.Errorf("failed to solve assembled linear system: %w", errSolve)
	}

	// Computes reaction forces [r_1] = (-1)·[k_11 d_1 + k_12 d_2] for Dirichlet-constrained dofs:
	scratch.Reset()
	scratch.ReuseAsVec(s.constrained)
	scratch.MulVec(k12, d2)
	r1.ScaleVec(-1.0, r1) // Turn Neumann loads into reactions
	r1.AddVec(r1, scratch)
	scratch.MulVec(k11, d1)
	r1.AddVec(r1, scratch)

	for _, transform := range p.EqTransforms {
		transform.Post(indices, r, d)
	}

	// Group primary/secondary solution with symbolic index
	for i := 0; i < s.dim; i++ {
		idx := indices.Unmap(i)
		s.d[i] = NodalValue{Index: idx, Value: d.AtVec(i)}
		s.r[i] = NodalValue{Index: idx, Value: r.AtVec(i)}
	}

	return s.d, s.r, indices.Failure()
}

func (s *linearSolver) initialise(dim, constrained int) error {
	if dim == constrained {
		return fmt.Errorf("all %v degrees of freedom have Dirichlet BC, no need to solve this", dim)
	}

	s.dim = dim
	s.constrained = constrained

	s.eqn = formMatrices(dim, constrained, &s.eqn)
	s.d = slices.Grow(s.d, dim)[0:dim]
	s.r = slices.Grow(s.r, dim)[0:dim]

	return nil
}

// formMatrices allocates or re-shapes existing matrices and defines the views over subsets of the
// backing matrices. When possible, existing buffers are re-used.
func formMatrices(dim, constrained int, prior *matrices) matrices {
	eqn := prior

	if eqn.k == nil {
		eqn.k = mat.NewSymDense(dim, nil)
		eqn.r = mat.NewVecDense(dim, nil)
		eqn.d = mat.NewVecDense(dim, nil)
		eqn.k12 = mat.NewDense(constrained, dim-constrained, nil)
		eqn.scratch = mat.NewVecDense(max(dim-constrained, constrained), nil)
	}

	for _, matrix := range []mat.Reseter{eqn.k, eqn.r, eqn.d, eqn.k12, eqn.scratch} {
		matrix.Reset()
	}

	eqn.k.ReuseAsSym(dim)
	eqn.r.ReuseAsVec(dim)
	eqn.d.ReuseAsVec(dim)
	eqn.k12.ReuseAs(constrained, dim-constrained)
	eqn.scratch.ReuseAsVec(max(dim-constrained, constrained))

	// The system k·d = r is partitioned as
	//   ⎡k_11 k_12⎤⎡d_1⎤  ⎡r_1⎤
	//   ⎣k_21 k_22⎦⎣d_2⎦  ⎣r_2⎦
	// where d_1 are all Dirichlet-prescribed values, and d_2 are free primary nodal values to be
	// solved for. We solve these through the linear system
	//   [k_22][d_2] = [r_2 - k_21 d_1]
	// and then use the solution vector to compute the reaction forces:
	//   [r_1] = (-1)·[k_11 d_1 + k_12 d_2]
	eqn.k11 = eqn.k.SliceSym(0, constrained)
	eqn.k22 = eqn.k.SliceSym(constrained, dim)

	eqn.d1 = eqn.d.SliceVec(0, constrained).(*mat.VecDense)
	eqn.d2 = eqn.d.SliceVec(constrained, dim).(*mat.VecDense)
	eqn.r1 = eqn.r.SliceVec(0, constrained).(*mat.VecDense)
	eqn.r2 = eqn.r.SliceVec(constrained, dim).(*mat.VecDense)

	return *eqn
}

func (s *linearSolver) extractK12() {
	for i := 0; i < s.constrained; i++ {
		for j := s.constrained; j < s.dim; j++ {
			s.eqn.k12.Set(i, j-s.constrained, s.eqn.k.At(i, j))
		}
	}
}
