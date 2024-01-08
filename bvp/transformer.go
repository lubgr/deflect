package bvp

import (
	"fmt"
	"log"
	"math"

	"github.com/lubgr/deflect/xyz"
	"gonum.org/v1/gonum/mat"
)

// Transformer mutate the system of equation before and after it is solved. For example, an inclined
// support can be split into a Transformer and a Dirichlet BC, such that they act independently to
// link degrees of freedom through a trigonometric relation. Transformer instances don't prescribe
// values in rhs or primary (this is done by NodalBC instances).
type Transformer interface {
	Pre(indices map[xyz.Index]int, k *mat.SymDense, rhs, primary *mat.VecDense) error
	Post(indices map[xyz.Index]int, rhs, primary *mat.VecDense) error
}

// NewAngularLinkTransformer instantiates a transformer to represent an inclined support to link
// displacements through an angle.
func NewAngularLinkTransformer(from, to xyz.Index, angle float64) Transformer {
	return &angularLink{from: from, to: to, c: math.Cos(angle), s: math.Sin(angle), irow: nil}
}

type angularLink struct {
	from, to   xyz.Index
	c, s       float64
	irow, jrow *mat.VecDense
}

type angularLinkPhase uint8

const (
	angularLinkPhasePre angularLinkPhase = iota
	angularLinkPhasePost
)

func (l *angularLink) Pre(
	indices map[xyz.Index]int,
	k *mat.SymDense,
	rhs, primary *mat.VecDense,
) error {

	// Initialise the buffer if necessary
	if l.irow == nil {
		dim := k.SymmetricDim()
		l.irow = mat.NewVecDense(dim, nil)
		l.jrow = mat.NewVecDense(dim, nil)
	}

	i, j, err := l.mapWithErr(indices)

	if err != nil {
		return fmt.Errorf("apply pre-step of equation transformation: %w", err)
	}

	l.transformTangent(i, j, k)
	l.transformVec(i, j, angularLinkPhasePre, rhs)
	l.transformVec(i, j, angularLinkPhasePre, primary)

	return nil
}

func (l *angularLink) Post(indices map[xyz.Index]int, rhs, primary *mat.VecDense) error {
	// No need to check for the initialisation of l's bufffers, since l.Pre must have been called
	// before l.Post, and l.Pre would initialise them.
	i, j, err := l.mapWithErr(indices)

	if err != nil {
		return fmt.Errorf("couldn't apply post-transformation: %w", err)
	}

	l.transformVec(i, j, angularLinkPhasePost, rhs)
	l.transformVec(i, j, angularLinkPhasePost, primary)

	return nil
}

func (l *angularLink) mapWithErr(indices map[xyz.Index]int) (from, to int, err error) {
	var ok bool

	if from, ok = indices[l.from]; !ok {
		return -1, -1, fmt.Errorf("index lookup for %v failed", l.from)
	}

	if to, ok = indices[l.to]; !ok {
		return -1, -1, fmt.Errorf("index lookup for %v failed", l.to)
	}

	return from, to, nil
}

func (l *angularLink) transformVec(i, j int, phase angularLinkPhase, v *mat.VecDense) {
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

func (l *angularLink) transformTangent(i, j int, k *mat.SymDense) {
	dim := k.SymmetricDim()
	s, c := l.s, l.c

	kii := k.At(i, i)
	kij := k.At(i, j)
	kjj := k.At(j, j)

	for m := range dim {
		if m == i || m == j {
			continue
		}
		l.irow.SetVec(m, c*k.At(m, i)+s*k.At(m, j))
		l.jrow.SetVec(m, -s*k.At(m, i)+c*k.At(m, j))
	}

	for m := range dim {
		k.SetSym(i, m, l.irow.AtVec(m))
		k.SetSym(j, m, l.jrow.AtVec(m))
	}

	k.SetSym(i, i, c*(c*kii+s*kij)+s*(c*kij+s*kjj))
	k.SetSym(i, j, c*(c*kij+s*kjj)-s*(c*kii+s*kij))
	k.SetSym(j, j, c*(c*kjj-s*kij)-s*(c*kij-s*kii))
}
