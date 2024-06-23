package deflect

import (
	"errors"

	"gonum.org/v1/gonum/mat"
)

// NewFrame3d returns a 3d beam element implementation.
func NewFrame3d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	common, errCommon := newOneDimElement(id, n0, n1, material)
	truss, errTruss := NewTruss3d(id, n0, n1, material, hinges)
	// TODO  newBeam3d
	beam, errBeam := newBeam2d(id, n0, n1, material, hinges)

	return &frame{
			oneDimElement: common,
			truss:         truss,
			beam:          beam,
		}, errors.Join(
			errCommon,
			errTruss,
			errBeam,
		)
}

// NewFrame2d returns a 2d beam element implementation.
func NewFrame2d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	common, errCommon := newOneDimElement(id, n0, n1, material)
	truss, errTruss := NewTruss2d(id, n0, n1, material, hinges)
	beam, errBeam := newBeam2d(id, n0, n1, material, hinges)

	return &frame{
			oneDimElement: common,
			truss:         truss,
			beam:          beam,
		}, errors.Join(
			errCommon,
			errTruss,
			errBeam,
		)
}

type frame struct {
	oneDimElement
	truss, beam Element
}

func (f *frame) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	f.truss.Assemble(indices, k, r, d)
	f.beam.Assemble(indices, k, r, d)
}

func (f *frame) Indices(set map[Index]struct{}) {
	f.truss.Indices(set)
	f.beam.Indices(set)
}

func (f *frame) AddLoad(bc NeumannElementBC) bool {
	return f.truss.AddLoad(bc) || f.beam.AddLoad(bc)
}

func (f *frame) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	s0 := f.truss.Interpolate(indices, which, d)
	s1 := f.beam.Interpolate(indices, which, d)

	return append(s0, s1...)
}
