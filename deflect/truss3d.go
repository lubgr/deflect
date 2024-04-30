package deflect

import (
	"errors"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type truss3d struct {
	// Much logic is shared between the 2d and the 3d implementation. We delegate to a 2d instance
	// where possible:
	truss2d
}

// NewTruss3d returns a new 3d truss implementation.
func NewTruss3d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	base, errTruss2d := NewTruss2d(id, n0, n1, material, hinges)

	if errTruss2d != nil {
		return nil, fmt.Errorf("failed to instantiate new 3d truss: %w", errors.Unwrap(errTruss2d))
	}

	concrete, ok := base.(*truss2d)

	if !ok {
		return nil, errors.New("bug: can't downcast fresh truss2d instance")
	}

	return &truss3d{truss2d: *concrete}, nil
}

func (t *truss3d) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	kAdd := func(i, j int, value float64) {
		k.SetSym(i, j, k.At(i, j)+value)
	}
	rAdd := func(i int, value float64) {
		r.SetVec(i, r.AtVec(i)+value)
	}

	idx := t.indicesAsArray()
	ux0, uy0, uz0 := indices.mapThree(idx[0], idx[1], idx[2])
	ux1, uy1, uz1 := indices.mapThree(idx[3], idx[4], idx[5])

	l := length(t.n0, t.n1)
	cx, cy, cz := directionCosine3d(t.n0, t.n1)

	kl := t.localNoHingeTangent(l)
	rl := t.localNoHingeLoads(l)

	t.hinges.reduce(kl, rl)

	kAdd(ux0, ux0, cx*cx*kl.At(0, 0))
	kAdd(ux0, uy0, cx*cy*kl.At(0, 0))
	kAdd(ux0, uz0, cx*cz*kl.At(0, 0))
	kAdd(ux0, ux1, cx*cx*kl.At(0, 1))
	kAdd(ux0, uy1, cx*cy*kl.At(0, 1))
	kAdd(ux0, uz1, cx*cz*kl.At(0, 1))

	kAdd(uy0, uy0, cy*cy*kl.At(0, 0))
	kAdd(uy0, uz0, cy*cz*kl.At(0, 0))
	kAdd(uy0, ux1, cx*cy*kl.At(0, 1))
	kAdd(uy0, uy1, cy*cy*kl.At(0, 1))
	kAdd(uy0, uz1, cy*cz*kl.At(0, 1))

	kAdd(uz0, uz0, cz*cz*kl.At(0, 0))
	kAdd(uz0, ux1, cx*cz*kl.At(0, 1))
	kAdd(uz0, uy1, cy*cz*kl.At(0, 1))
	kAdd(uz0, uz1, cz*cz*kl.At(0, 1))

	kAdd(ux1, ux1, cx*cx*kl.At(1, 1))
	kAdd(ux1, uy1, cx*cy*kl.At(1, 1))
	kAdd(ux1, uz1, cx*cz*kl.At(1, 1))

	kAdd(uy1, uy1, cy*cy*kl.At(1, 1))
	kAdd(uy1, uz1, cy*cz*kl.At(1, 1))

	kAdd(uz1, uz1, cz*cz*kl.At(1, 1))

	rAdd(ux0, cx*rl.AtVec(0))
	rAdd(uy0, cy*rl.AtVec(0))
	rAdd(uz0, cz*rl.AtVec(0))

	rAdd(ux1, cx*rl.AtVec(1))
	rAdd(uy1, cy*rl.AtVec(1))
	rAdd(uz1, cz*rl.AtVec(1))
}

func (t *truss3d) indicesAsArray() *[6]Index {
	indices := [...]Index{
		{NodalID: t.n0.ID, Dof: Ux},
		{NodalID: t.n0.ID, Dof: Uy},
		{NodalID: t.n0.ID, Dof: Uz},
		{NodalID: t.n1.ID, Dof: Ux},
		{NodalID: t.n1.ID, Dof: Uy},
		{NodalID: t.n1.ID, Dof: Uz},
	}

	return &indices
}

func (t *truss3d) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	switch which {
	case FctUx, FctNx:
	default:
		return nil
	}

	ux0, nx0 := t.startNodeValues(indices, d)
	return t.InterpolateNxOrUx(which, ux0, nx0)
}

func (t *truss3d) startNodeValues(indices EqLayout, d *mat.VecDense) (dx0, nx0 float64) {
	l := length(t.n0, t.n1)
	cx, cy, cz := directionCosine3d(t.n0, t.n1)

	kl := t.localNoHingeTangent(l)
	rl := t.localNoHingeLoads(l)

	idx := t.indicesAsArray()
	ux0, uy0, uz0 := indices.mapThree(idx[0], idx[1], idx[2])
	ux1, uy1, uz1 := indices.mapThree(idx[3], idx[4], idx[5])

	d0, d1, d2 := d.AtVec(ux0), d.AtVec(uy0), d.AtVec(uz0)
	d3, d4, d5 := d.AtVec(ux1), d.AtVec(uy1), d.AtVec(uz1)

	dl := mat.NewVecDense(2, []float64{
		cx*d0 + cy*d1 + cz*d2,
		cx*d3 + cy*d4 + cz*d5,
	})

	t.hinges.enhance(kl, rl, dl)

	dx0 = dl.AtVec(0)

	dl.MulVec(kl, dl) // Stores local end forces/stresses now
	rl.SubVec(rl, dl)

	nx0 = rl.AtVec(0)

	return dx0, nx0
}

func (t *truss3d) Indices(set map[Index]struct{}) {
	for _, index := range t.indicesAsArray() {
		set[index] = struct{}{}
	}
}
