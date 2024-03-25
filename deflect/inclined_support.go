package deflect

import (
	"log"
	"math"

	"gonum.org/v1/gonum/mat"
)

// NewInclinedSupport instantiates a transformer and a zero Dirichlet boundary condition. Together,
// these two objects represent an inclined support. They are tightly coupled, as the prescribed
// index determines how the Transformer acts on the assembled tangent and vectors, and hence
// returned by the same constructor(-like) function.
func NewInclinedSupport(from, to Index, angle float64) (Transformer, NodalValue) {
	return &inclinedSupport{from: from, to: to, c: math.Cos(angle), s: math.Sin(angle), irow: nil},
		NodalValue{Index: to, Value: 0}
}

type inclinedSupport struct {
	from, to   Index
	c, s       float64
	irow, jrow *mat.VecDense
}

type angularLinkPhase uint8

const (
	angularLinkPhasePre angularLinkPhase = iota
	angularLinkPhasePost
)

func (l *inclinedSupport) Pre(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	// Initialise the buffer if necessary. We might want to think about sharing the same buffer
	// between multiple inclinedSupport instances at some point, since re-using the buffer would be
	// more efficient.
	if l.irow == nil {
		dim := k.SymmetricDim()
		l.irow = mat.NewVecDense(dim, nil)
		l.jrow = mat.NewVecDense(dim, nil)
	}

	i, j := indices.mapTwo(l.from, l.to)

	l.transformTangent(i, j, k)
	l.transformVec(i, j, angularLinkPhasePre, r)
	l.transformVec(i, j, angularLinkPhasePre, d)
}

func (l *inclinedSupport) Post(indices EqLayout, r, d *mat.VecDense) {
	// No need to check for the initialisation of l's bufffers, since l.Pre must have been called
	// before l.Post, and l.Pre would initialise them.
	i, j := indices.mapTwo(l.from, l.to)

	l.transformVec(i, j, angularLinkPhasePost, r)
	l.transformVec(i, j, angularLinkPhasePost, d)
}

func (l *inclinedSupport) transformVec(i, j int, phase angularLinkPhase, v *mat.VecDense) {
	// Computes tᵀ·v for the Pre-step, and t*v for the Post-step, where v is the given vector, and t
	// is the coordinate transformation
	//   ⎡c  -s⎤
	//   ⎣s   c⎦
	s, c := l.s, l.c
	vi := v.AtVec(i)
	vj := v.AtVec(j)

	switch phase {
	case angularLinkPhasePre:
		v.SetVec(i, c*vi+s*vj)
		v.SetVec(j, -s*vi+c*vj)
	case angularLinkPhasePost:
		v.SetVec(i, c*vi-s*vj)
		v.SetVec(j, s*vi+c*vj)
	default:
		log.Printf("Unknown enumerator %v for angular link phase", phase)
	}
}

func (l *inclinedSupport) transformTangent(i, j int, k *mat.SymDense) {
	// Computes tᵀ·k·t where t is the coordinate transformation
	//   ⎡c  -s⎤
	//   ⎣s   c⎦
	// This implementation tries to use as little storage as possible. Arguably, this could also be
	// done by using using a sparse matrix library. But this is left as an optimisation exercise (that
	// should demonstrate the computation benefit over this hand-rolled scheme).
	dim := k.SymmetricDim()
	s, c := l.s, l.c

	kii := k.At(i, i)
	kij := k.At(i, j)
	kjj := k.At(j, j)

	for m := range dim {
		if m == i || m == j {
			continue
		}
		// Populate scratch buffers for the matrix-matrix multiplication that we carry out only for the
		// two relevant rows.
		l.irow.SetVec(m, c*k.At(m, i)+s*k.At(m, j))
		l.jrow.SetVec(m, -s*k.At(m, i)+c*k.At(m, j))
	}

	// Spill scratch buffers into the destination matrix
	for m := range dim {
		k.SetSym(i, m, l.irow.AtVec(m))
		k.SetSym(j, m, l.jrow.AtVec(m))
	}

	// Finally, the diagonal terms
	k.SetSym(i, i, c*(c*kii+s*kij)+s*(c*kij+s*kjj))
	k.SetSym(i, j, c*(c*kij+s*kjj)-s*(c*kii+s*kij))
	k.SetSym(j, j, c*(c*kjj-s*kij)-s*(c*kij-s*kii))
}
