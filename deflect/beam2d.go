package deflect

import (
	"errors"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

type beam2d struct {
	oneDimElement
	hinges condenser
}

// newBeam2d returns a 2d beam element implementation.
func newBeam2d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	localIndices := [...]Index{
		{NodalID: n0.ID, Dof: Uz},
		{NodalID: n0.ID, Dof: Phiy},
		{NodalID: n1.ID, Dof: Uz},
		{NodalID: n1.ID, Dof: Phiy},
	}

	common, errCommon := newOneDimElement(id, n0, n1, material)
	condenser, errCondensation := newStaticHingeCondenser(localIndices[:], hinges)

	if err := errors.Join(errCommon, errCondensation); err != nil {
		return nil, fmt.Errorf("failed to instantiate new 2d beam: %w", err)
	}

	return &beam2d{oneDimElement: common, hinges: condenser}, nil
}

func (b *beam2d) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	kAdd := func(i, j int, value float64) {
		k.SetSym(i, j, k.At(i, j)+value)
	}
	rAdd := func(i int, value float64) {
		r.SetVec(i, r.AtVec(i)+value)
	}

	idx := b.indicesAsArray()
	ux0, uz0, phiy0 := indices.mapThree(idx[0], idx[1], idx[2])
	ux1, uz1, phiy1 := indices.mapThree(idx[3], idx[4], idx[5])

	s, c := sineCosine2d(b.n0, b.n1)
	s2, c2, sc := s*s, c*c, s*c

	l := length(b.n0, b.n1)
	kl := b.localNoHingeTangent(l)
	rl := b.localNoHingeLoads(l)

	b.hinges.reduce(kl, rl)

	kAdd(ux0, ux0, s2*kl.At(0, 0))
	kAdd(ux0, uz0, -sc*kl.At(0, 0))
	kAdd(ux0, phiy0, -s*kl.At(0, 1))
	kAdd(ux0, ux1, s2*kl.At(0, 2))
	kAdd(ux0, uz1, -sc*kl.At(0, 2))
	kAdd(ux0, phiy1, -s*kl.At(0, 3))

	kAdd(uz0, uz0, c2*kl.At(0, 0))
	kAdd(uz0, phiy0, c*kl.At(0, 1))
	kAdd(uz0, ux1, -sc*kl.At(0, 2))
	kAdd(uz0, uz1, c2*kl.At(0, 2))
	kAdd(uz0, phiy1, c*kl.At(0, 3))

	kAdd(phiy0, phiy0, kl.At(1, 1))
	kAdd(phiy0, ux1, -s*kl.At(1, 2))
	kAdd(phiy0, uz1, c*kl.At(1, 2))
	kAdd(phiy0, phiy1, kl.At(1, 3))

	kAdd(ux1, ux1, s2*kl.At(2, 2))
	kAdd(ux1, uz1, -sc*kl.At(2, 2))
	kAdd(ux1, phiy1, -s*kl.At(2, 3))

	kAdd(uz1, uz1, c2*kl.At(2, 2))
	kAdd(uz1, phiy1, c*kl.At(2, 3))

	kAdd(phiy1, phiy1, kl.At(3, 3))

	rAdd(ux0, s*rl.AtVec(0))
	rAdd(uz0, -c*rl.AtVec(0))
	rAdd(phiy0, -rl.AtVec(1))
	rAdd(ux1, s*rl.AtVec(2))
	rAdd(uz1, -c*rl.AtVec(2))
	rAdd(phiy1, -rl.AtVec(3))
}

func (b *beam2d) localNoHingeTangent(l float64) *mat.SymDense {
	k := mat.NewSymDense(4, nil)
	EI := b.material.YoungsModulus * b.material.Iyy()
	l2, l3 := l*l, l*l*l

	k.SetSym(0, 0, 12*EI/l3)
	k.SetSym(0, 1, -6*EI/l2)
	k.SetSym(0, 2, -12*EI/l3)
	k.SetSym(0, 3, -6*EI/l2)

	k.SetSym(1, 1, 4*EI/l)
	k.SetSym(1, 2, 6*EI/l2)
	k.SetSym(1, 3, 2*EI/l)

	k.SetSym(2, 2, 12*EI/l3)
	k.SetSym(2, 3, 6*EI/l2)

	k.SetSym(3, 3, 4*EI/l)

	return k
}

func (b *beam2d) localNoHingeLoads(l float64) *mat.VecDense {
	var rz0, rphiy0, rz1, rphiy1 float64

	for _, bc := range b.loads {
		loadDispatch(bc,
			func(load *neumannConcentrated) {
				if load.kind == Uz {
					a, b, fz := load.position/l, 1.0-load.position/l, load.value
					rz0 += fz * b * b * (3 - 2*b)
					rz1 += fz * a * a * (3 - 2*a)
					rphiy0 -= fz * a * b * b * l
					rphiy1 += fz * b * a * a * l
				} else if load.kind == Phiy {
					a, b, my := load.position/l, 1.0-load.position/l, load.value
					rz0 += my * 6 / l * a * b
					rz1 -= my * 6 / l * a * b
					rphiy0 -= my * b * (3*a - 1)
					rphiy1 -= my * a * (3*b - 1)
				}
			},
			func(load *neumannConstant) {
				if load.kind == Uz {
					q := load.value
					rz0 += q * l / 2
					rz1 += q * l / 2
					rphiy0 -= q * l * l / 12
					rphiy1 += q * l * l / 12
				}
			},
			func(load *neumannLinear) {
				if load.kind == Uz {
					q0, q1 := load.first, load.last
					rz0 += (7*q0 + 3*q1) * l / 20
					rz1 += (3*q0 + 7*q1) * l / 20
					rphiy0 -= (3*q0 + 2*q1) * l * l / 60
					rphiy1 += (2*q0 + 3*q1) * l * l / 60
				}
			})
	}

	r := mat.NewVecDense(4, nil)
	r.SetVec(0, rz0)
	r.SetVec(1, rphiy0)
	r.SetVec(2, rz1)
	r.SetVec(3, rphiy1)

	return r
}

func (b *beam2d) indicesAsArray() *[6]Index {
	indices := [...]Index{
		{NodalID: b.n0.ID, Dof: Ux},
		{NodalID: b.n0.ID, Dof: Uz},
		{NodalID: b.n0.ID, Dof: Phiy},
		{NodalID: b.n1.ID, Dof: Ux},
		{NodalID: b.n1.ID, Dof: Uz},
		{NodalID: b.n1.ID, Dof: Phiy},
	}

	return &indices
}

func (b *beam2d) AddLoad(bc NeumannElementBC) bool {
	supported := false

	loadDispatch(bc,
		func(l *neumannConcentrated) { supported = l.kind == Uz || l.kind == Phiy },
		func(l *neumannConstant) { supported = l.kind == Uz },
		func(l *neumannLinear) { supported = l.kind == Uz })

	return supported && b.oneDimElement.AddLoad(bc)
}

func (b *beam2d) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	switch which {
	case FctVz, FctMy, FctPhiy, FctUz:
	default:
		return nil
	}

	// All interpolations are based on My; Vz by differentiation, phiy and uz by integration. Vz,
	// phiy, and uz don't have to be computed like that, but it seemed like a good strategy; a single
	// interpolation, My, must be understood and implemented for every possible element loading. Then,
	// all other interpolations are derived from that.
	my := b.InterpolateMy(indices, d)

	switch which {
	case FctMy:
		return my
	case FctVz:
		return transform(func(p PolyPiece) PolyPiece { return p.derive() }, my)
	}

	// In what follows

	EI := b.material.YoungsModulus * b.material.Iyy()
	uz0, phiy0, _, _ := b.startNodeValues(indices, d)
	myOverEI := my // Just a renaming, coefficients are shallow-copied
	myOverEI.multiply(1 / EI)

	phiy := my.integrate(phiy0)
	// Given the local co-system, we have d/dx w(x) = -phiy(x). We integrate first without this sign
	// flip, but since uz0 is already in the local co-system, we need to invert its sign and then
	// invert the sign of the entire polynomial altogether afterwards.
	uz := phiy.integrate(-uz0)
	uz.multiply(-1)

	switch which {
	case FctPhiy:
		return phiy
	case FctUz:
		return uz
	}

	return nil
}

func (b *beam2d) InterpolateMy(indices EqLayout, d *mat.VecDense) PolySequence {
	_, _, my0, vz0 := b.startNodeValues(indices, d)
	l := length(b.n0, b.n1)
	result := PolySequence{PolyPiece{X0: 0, XE: l, Coeff: []float64{my0, vz0}}}

	for _, bc := range b.loads {
		var p PolyPiece

		loadDispatch(bc,
			func(load *neumannConcentrated) {
				if load.kind == Uz {
					fz, a := load.value, load.position
					p = PolyPiece{X0: a, XE: l, Coeff: []float64{fz * a, -fz}}
				} else if load.kind == Phiy {
					my, a := load.value, load.position
					p = PolyPiece{X0: a, XE: l, Coeff: []float64{-my}}
				}
			},
			func(load *neumannConstant) {
				if load.kind == Uz {
					q := load.value
					p = PolyPiece{X0: 0, XE: l, Coeff: []float64{0, 0, -q / 2}}
				}
			},
			func(load *neumannLinear) {
				if load.kind == Uz {
					q0, qE := load.first, load.last
					p = PolyPiece{X0: 0, XE: l, Coeff: []float64{0, 0, -0.5 * q0, -(qE - q0) / (6 * l)}}
				}
			})

		result = append(result, p)
	}

	return result.flatten()
}

func (b *beam2d) startNodeValues(indices EqLayout, d *mat.VecDense) (uz, phiy, my, vz float64) {
	l := length(b.n0, b.n1)
	s, c := sineCosine2d(b.n0, b.n1)

	kl := b.localNoHingeTangent(l)
	rl := b.localNoHingeLoads(l)

	idx := b.indicesAsArray()
	ux0, uz0, phiy0 := indices.mapThree(idx[0], idx[1], idx[2])
	ux1, uz1, phiy1 := indices.mapThree(idx[3], idx[4], idx[5])

	d0, d1, d2 := d.AtVec(ux0), d.AtVec(uz0), d.AtVec(phiy0)
	d3, d4, d5 := d.AtVec(ux1), d.AtVec(uz1), d.AtVec(phiy1)
	dl := mat.NewVecDense(4, []float64{
		s*d0 - c*d1,
		-d2,
		s*d3 - c*d4,
		-d5,
	})

	b.hinges.enhance(kl, rl, dl)

	uz, phiy = dl.AtVec(0), dl.AtVec(1)

	dl.MulVec(kl, dl) // Stores local end forces/stresses now
	rl.SubVec(rl, dl)

	vz, my = rl.AtVec(0), rl.AtVec(1)

	return
}

func (b *beam2d) Indices(set map[Index]struct{}) {
	for _, index := range b.indicesAsArray() {
		set[index] = struct{}{}
	}
}
