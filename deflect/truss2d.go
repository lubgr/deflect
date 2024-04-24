package deflect

import (
	"errors"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type truss2d struct {
	oneDimElement
	hinges condenser
}

// NewTruss2d returns a new linear-elastic 2d truss implementation.
func NewTruss2d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	localIndices := [...]Index{{NodalID: n0.ID, Dof: Ux}, {NodalID: n1.ID, Dof: Ux}}

	common, errCommon := newOneDimElement(id, n0, n1, material)
	condenser, errCondensation := newStaticHingeCondenser(localIndices[:], hinges)

	if err := errors.Join(errCommon, errCondensation); err != nil {
		return nil, fmt.Errorf("failed to instantiate new 2d truss: %w", err)
	}

	return &truss2d{oneDimElement: common, hinges: condenser}, nil
}

func (t *truss2d) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	kAdd := func(i, j int, value float64) {
		k.SetSym(i, j, k.At(i, j)+value)
	}
	rAdd := func(i int, value float64) {
		r.SetVec(i, r.AtVec(i)+value)
	}

	idx := t.indicesAsArray()
	ux0, uz0, ux1, uz1 := indices.mapFour(idx[0], idx[1], idx[2], idx[3])

	s, c := sineCosine2d(t.n0, t.n1)
	s2, c2, sc := s*s, c*c, s*c
	l := length(t.n0, t.n1)
	kl := t.localNoHingeTangent(l)
	rl := t.localNoHingeLoads(l)

	t.hinges.reduce(kl, rl)

	kAdd(ux0, ux0, c2*kl.At(0, 0))
	kAdd(ux0, uz0, sc*kl.At(0, 0))
	kAdd(ux0, ux1, c2*kl.At(0, 1))
	kAdd(ux0, uz1, sc*kl.At(0, 1))

	kAdd(uz0, uz0, s2*kl.At(0, 0))
	kAdd(uz0, ux1, sc*kl.At(0, 1))
	kAdd(uz0, uz1, s2*kl.At(0, 1))

	kAdd(ux1, ux1, c2*kl.At(1, 1))
	kAdd(ux1, uz1, sc*kl.At(1, 1))

	kAdd(uz1, uz1, s2*kl.At(1, 1))

	rAdd(ux0, c*rl.AtVec(0))
	rAdd(uz0, s*rl.AtVec(0))
	rAdd(ux1, c*rl.AtVec(1))
	rAdd(uz1, s*rl.AtVec(1))
}

func (t *truss2d) localNoHingeTangent(l float64) *mat.SymDense {
	k := mat.NewSymDense(2, nil)
	EA := t.material.YoungsModulus * t.material.Area()

	k.SetSym(0, 0, EA/l)
	k.SetSym(0, 1, -EA/l)
	k.SetSym(1, 1, EA/l)

	return k
}

func (t *truss2d) localNoHingeLoads(l float64) *mat.VecDense {
	var rx0, rx1 float64

	for _, bc := range t.loads {
		loadDispatch(bc,
			func(load *neumannElementConcentrated) {
				a, b, fx := load.position, l-load.position, load.value
				rx0 += fx * b / l
				rx1 += fx * a / l
			},
			func(load *neumannElementConstant) {
				q := load.value
				rx0 += q * l / 2
				rx1 += q * l / 2
			},
			func(load *neumannElementLinear) {
				q0, q1 := load.first, load.last
				rx0 += l * (2*q0 + q1) / 6.0
				rx1 += l * (q0 + 2*q1) / 6.0
			})
	}

	r := mat.NewVecDense(2, nil)
	r.SetVec(0, rx0)
	r.SetVec(1, rx1)

	return r
}

func (t *truss2d) indicesAsArray() *[4]Index {
	indices := [...]Index{
		{NodalID: t.n0.ID, Dof: Ux},
		{NodalID: t.n0.ID, Dof: Uz},
		{NodalID: t.n1.ID, Dof: Ux},
		{NodalID: t.n1.ID, Dof: Uz},
	}

	return &indices
}

func (t *truss2d) Indices(set map[Index]struct{}) {
	for _, index := range t.indicesAsArray() {
		set[index] = struct{}{}
	}
}

func (t *truss2d) AddLoad(bc NeumannElementBC) bool {
	var kind Dof

	loadDispatch(bc,
		func(load *neumannElementConcentrated) {
			kind = load.kind
		},
		func(load *neumannElementConstant) {
			kind = load.kind
		},
		func(load *neumannElementLinear) {
			kind = load.kind
		})

	if kind == Ux {
		return t.oneDimElement.AddLoad(bc)
	}

	return false
}

func (t *truss2d) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	if which == FctNx {
		return t.InterpolateNx(indices, d)
	} else if which != FctUx {
		return nil
	}

	// The horizontal displacement is uₓ = ∫ εₓₓ dx = ∫ Nₓ/EA dx. We could also compute this
	// differently and dispatch over the element loads using manually integrated formulae. This
	// solution without another load-specific approach seems preferable for now.
	EA := t.material.YoungsModulus * t.material.Area()
	ux0, _ := t.startNodeValues(indices, d)

	eps := t.InterpolateNx(indices, d)
	eps.multiply(1 / EA)
	ux := eps.integrate(ux0)

	return ux
}

func (t *truss2d) InterpolateNx(indices EqLayout, d *mat.VecDense) PolySequence {
	_, nx0 := t.startNodeValues(indices, d)
	l := length(t.n0, t.n1)
	result := PolySequence{}

	result = append(result, PolyPiece{X0: 0, XE: l, Coeff: []float64{nx0}})

	for _, bc := range t.loads {
		// We accepted only Ux loads in AddLoad, so we don't have to inspect the load's kind again
		loadDispatch(bc,
			func(load *neumannElementConcentrated) {
				fx, a := load.value, load.position
				result = append(result, PolyPiece{X0: a, XE: l, Coeff: []float64{-fx}})
			},
			func(load *neumannElementConstant) {
				qx := load.value
				result = append(result, PolyPiece{X0: 0, XE: l, Coeff: []float64{0, -qx}})
			},
			func(load *neumannElementLinear) {
				q0, qE := load.first, load.last
				result = append(
					result,
					PolyPiece{X0: 0, XE: l, Coeff: []float64{0, -q0, -(qE - q0) / (2 * l)}},
				)
			})
	}

	return result.flatten()
}

func (t *truss2d) startNodeValues(indices EqLayout, d *mat.VecDense) (dx0, nx0 float64) {
	l := length(t.n0, t.n1)
	s, c := sineCosine2d(t.n0, t.n1)

	kl := t.localNoHingeTangent(l)
	rl := t.localNoHingeLoads(l)

	idx := t.indicesAsArray()
	ux0, uz0, ux1, uz1 := indices.mapFour(idx[0], idx[1], idx[2], idx[3])

	d0, d1, d2, d3 := d.AtVec(ux0), d.AtVec(uz0), d.AtVec(ux1), d.AtVec(uz1)
	dl := mat.NewVecDense(2, []float64{
		c*d0 + s*d1,
		c*d2 + s*d3,
	})

	dx0 = dl.AtVec(0)

	t.hinges.enhance(kl, rl, dl)

	dl.MulVec(kl, dl) // Stores local end forces/stresses now
	rl.SubVec(rl, dl)

	nx0 = rl.AtVec(0)

	return dx0, nx0
}
