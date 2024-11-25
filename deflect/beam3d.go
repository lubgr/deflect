package deflect

import (
	"errors"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// newBeam3d returns a 3d beam element implementation.
func newBeam3d(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	localIndicesV := [...]Index{
		{NodalID: n0.ID, Dof: Uz},
		{NodalID: n0.ID, Dof: Phiy},
		{NodalID: n1.ID, Dof: Uz},
		{NodalID: n1.ID, Dof: Phiy},
	}
	localIndicesH := [...]Index{
		{NodalID: n0.ID, Dof: Uy},
		{NodalID: n0.ID, Dof: Phiz},
		{NodalID: n1.ID, Dof: Uy},
		{NodalID: n1.ID, Dof: Phiz},
	}

	common, errCommon := newOneDimElement(id, n0, n1, material)
	condenserV, errCondensationV := newStaticHingeCondenser(localIndicesV[:], hinges)
	condenserH, errCondensationH := newStaticHingeCondenser(localIndicesH[:], hinges)

	if err := errors.Join(errCommon, errCondensationV, errCondensationH); err != nil {
		return nil, fmt.Errorf("failed to instantiate new 3d beam: %w", err)
	}

	return &beam3d{
		oneDimElement: common,
		delegate:      beam2d{oneDimElement: common, hinges: nil},
		hingesV:       condenserV,
		hingesH:       condenserH,
	}, nil
}

type beam3d struct {
	oneDimElement
	delegate         beam2d
	hingesV, hingesH condenser
}

func (b *beam3d) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	kAdd := func(i, j int, value float64) {
		k.SetSym(i, j, k.At(i, j)+value)
	}
	rAdd := func(i int, value float64) {
		r.SetVec(i, r.AtVec(i)+value)
	}

	idx := b.indicesAsArray()
	ux0, uy0, uz0 := indices.mapThree(idx[0], idx[1], idx[2])
	ux1, uy1, uz1 := indices.mapThree(idx[6], idx[7], idx[8])
	phix0, phiy0, phiz0 := indices.mapThree(idx[3], idx[4], idx[5])
	phix1, phiy1, phiz1 := indices.mapThree(idx[9], idx[10], idx[11])

	l := length(b.n0, b.n1)
	klv := b.delegate.localNoHingeTangent(l, b.material.Iyy())
	klh := b.delegate.localNoHingeTangent(l, b.material.Izz())
	rlv := b.delegate.localNoHingeLoads(l, Uz, Phiy)
	rlh := b.delegate.localNoHingeLoads(l, Uy, Phiz)

	b.hingesV.reduce(klv, rlv)
	b.hingesH.reduce(klh, rlh)

	// FIXME we currently have a material instance in oneDimElement and one in the oneDimElement
	// instance embedded in `delegate`.
	// - Should the beam2d implementation already use the roll angle?! Do we support rolled cross
	//   sections for 2d frames? (why not?!) If yes, how does this influence the construction of the
	//   3d rotation matrix here?
	// - If roll angle in 2d is supported, couldn't that be implemented by transforming Iyy
	//   accordingly? We'd still have no change of primary y/z axes, only the CS constants change.
	//   But what is CrossSection.Iyy() supposed to do then? In 2d, return a roll-angle adjusted
	//   value, but in 3d, return it as is?
	//   => Possible solution: we know upfront from the context if this is a 2d or a 3d problem. If
	//   it's 2d, construct the CS instances such that the roll angle is zero but the instance knows
	//   that it's being rotated (e.g. for the Generic/'constants' CS, pre-compute rotated values). If
	//   it's 3d, the roll angle is left as is, as well as the constants.
	//   Note that this might still not be ideal since it uses global knowledge to influence local
	//   behaviour. Is there a different solution? How does other software solve this?
	_, _, _, t10, t11, t12, t20, t21, t22 := rotation3d(b.n0, b.n1, 0)

	kAdd(ux0, ux0, klh.At(0, 0)*t10*t10+klv.At(0, 0)*t20*t20)
	kAdd(ux0, uy0, klh.At(0, 0)*t10*t11+klv.At(0, 0)*t20*t21)
	kAdd(ux0, uz0, klh.At(0, 0)*t10*t12+klv.At(0, 0)*t20*t22)
	kAdd(ux0, phix0, klh.At(0, 1)*t10*t20+klv.At(0, 1)*t10*t20)
	kAdd(ux0, phiy0, klh.At(0, 1)*t10*t21+klv.At(0, 1)*t11*t20)
	kAdd(ux0, phiz0, klh.At(0, 1)*t10*t22+klv.At(0, 1)*t12*t20)
	kAdd(ux0, ux1, klh.At(0, 2)*t10*t10+klv.At(0, 2)*t20*t20)
	kAdd(ux0, uy1, klh.At(0, 2)*t10*t11+klv.At(0, 2)*t20*t21)
	kAdd(ux0, uz1, klh.At(0, 2)*t10*t12+klv.At(0, 2)*t20*t22)
	kAdd(ux0, phix1, klh.At(0, 3)*t10*t20+klv.At(0, 3)*t10*t20)
	kAdd(ux0, phiy1, klh.At(0, 3)*t10*t21+klv.At(0, 3)*t11*t20)
	kAdd(ux0, phiz1, klh.At(0, 3)*t10*t22+klv.At(0, 3)*t12*t20)

	kAdd(uy0, uy0, klh.At(0, 0)*t11*t11+klv.At(0, 0)*t21*t21)
	kAdd(uy0, uz0, klh.At(0, 0)*t11*t12+klv.At(0, 0)*t21*t22)
	kAdd(uy0, phix0, klh.At(0, 1)*t11*t20+klv.At(0, 1)*t10*t21)
	kAdd(uy0, phiy0, klh.At(0, 1)*t11*t21+klv.At(0, 1)*t11*t21)
	kAdd(uy0, phiz0, klh.At(0, 1)*t11*t22+klv.At(0, 1)*t12*t21)
	kAdd(uy0, ux1, klh.At(0, 2)*t10*t11+klv.At(0, 2)*t20*t21)
	kAdd(uy0, uy1, klh.At(0, 2)*t11*t11+klv.At(0, 2)*t21*t21)
	kAdd(uy0, uz1, klh.At(0, 2)*t11*t12+klv.At(0, 2)*t21*t22)
	kAdd(uy0, phix1, klh.At(0, 3)*t11*t20+klv.At(0, 3)*t10*t21)
	kAdd(uy0, phiy1, klh.At(0, 3)*t11*t21+klv.At(0, 3)*t11*t21)
	kAdd(uy0, phiz1, klh.At(0, 3)*t11*t22+klv.At(0, 3)*t12*t21)

	kAdd(uz0, uz0, klh.At(0, 0)*t12*t12+klv.At(0, 0)*t22*t22)
	kAdd(uz0, phix0, klh.At(0, 1)*t12*t20+klv.At(0, 1)*t10*t22)
	kAdd(uz0, phiy0, klh.At(0, 1)*t12*t21+klv.At(0, 1)*t11*t22)
	kAdd(uz0, phiz0, klh.At(0, 1)*t12*t22+klv.At(0, 1)*t12*t22)
	kAdd(uz0, ux1, klh.At(0, 2)*t10*t12+klv.At(0, 2)*t20*t22)
	kAdd(uz0, uy1, klh.At(0, 2)*t11*t12+klv.At(0, 2)*t21*t22)
	kAdd(uz0, uz1, klh.At(0, 2)*t12*t12+klv.At(0, 2)*t22*t22)
	kAdd(uz0, phix1, klh.At(0, 3)*t12*t20+klv.At(0, 3)*t10*t22)
	kAdd(uz0, phiy1, klh.At(0, 3)*t12*t21+klv.At(0, 3)*t11*t22)
	kAdd(uz0, phiz1, klh.At(0, 3)*t12*t22+klv.At(0, 3)*t12*t22)

	kAdd(phix0, phix0, klh.At(1, 1)*t20*t20+klv.At(1, 1)*t10*t10)
	kAdd(phix0, phiy0, klh.At(1, 1)*t20*t21+klv.At(1, 1)*t10*t11)
	kAdd(phix0, phiz0, klh.At(1, 1)*t20*t22+klv.At(1, 1)*t10*t12)
	kAdd(phix0, ux1, klh.At(1, 2)*t10*t20+klv.At(1, 2)*t10*t20)
	kAdd(phix0, uy1, klh.At(1, 2)*t11*t20+klv.At(1, 2)*t10*t21)
	kAdd(phix0, uz1, klh.At(1, 2)*t12*t20+klv.At(1, 2)*t10*t22)
	kAdd(phix0, phix1, klh.At(1, 3)*t20*t20+klv.At(1, 3)*t10*t10)
	kAdd(phix0, phiy1, klh.At(1, 3)*t20*t21+klv.At(1, 3)*t10*t11)
	kAdd(phix0, phiz1, klh.At(1, 3)*t20*t22+klv.At(1, 3)*t10*t12)

	kAdd(phiy0, phiy0, klh.At(1, 1)*t21*t21+klv.At(1, 1)*t11*t11)
	kAdd(phiy0, phiz0, klh.At(1, 1)*t21*t22+klv.At(1, 1)*t11*t12)
	kAdd(phiy0, ux1, klh.At(1, 2)*t10*t21+klv.At(1, 2)*t11*t20)
	kAdd(phiy0, uy1, klh.At(1, 2)*t11*t21+klv.At(1, 2)*t11*t21)
	kAdd(phiy0, uz1, klh.At(1, 2)*t12*t21+klv.At(1, 2)*t11*t22)
	kAdd(phiy0, phix1, klh.At(1, 3)*t20*t21+klv.At(1, 3)*t10*t11)
	kAdd(phiy0, phiy1, klh.At(1, 3)*t21*t21+klv.At(1, 3)*t11*t11)
	kAdd(phiy0, phiz1, klh.At(1, 3)*t21*t22+klv.At(1, 3)*t11*t12)

	kAdd(phiz0, phiz0, klh.At(1, 1)*t22*t22+klv.At(1, 1)*t12*t12)
	kAdd(phiz0, ux1, klh.At(1, 2)*t10*t22+klv.At(1, 2)*t12*t20)
	kAdd(phiz0, uy1, klh.At(1, 2)*t11*t22+klv.At(1, 2)*t12*t21)
	kAdd(phiz0, uz1, klh.At(1, 2)*t12*t22+klv.At(1, 2)*t12*t22)
	kAdd(phiz0, phix1, klh.At(1, 3)*t20*t22+klv.At(1, 3)*t10*t12)
	kAdd(phiz0, phiy1, klh.At(1, 3)*t21*t22+klv.At(1, 3)*t11*t12)
	kAdd(phiz0, phiz1, klh.At(1, 3)*t22*t22+klv.At(1, 3)*t12*t12)

	kAdd(ux1, ux1, klh.At(2, 2)*t10*t10+klv.At(2, 2)*t20*t20)
	kAdd(ux1, uy1, klh.At(2, 2)*t10*t11+klv.At(2, 2)*t20*t21)
	kAdd(ux1, uz1, klh.At(2, 2)*t10*t12+klv.At(2, 2)*t20*t22)
	kAdd(ux1, phix1, klh.At(2, 3)*t10*t20+klv.At(2, 3)*t10*t20)
	kAdd(ux1, phiy1, klh.At(2, 3)*t10*t21+klv.At(2, 3)*t11*t20)
	kAdd(ux1, phiz1, klh.At(2, 3)*t10*t22+klv.At(2, 3)*t12*t20)

	kAdd(uy1, uy1, klh.At(2, 2)*t11*t11+klv.At(2, 2)*t21*t21)
	kAdd(uy1, uz1, klh.At(2, 2)*t11*t12+klv.At(2, 2)*t21*t22)
	kAdd(uy1, phix1, klh.At(2, 3)*t11*t20+klv.At(2, 3)*t10*t21)
	kAdd(uy1, phiy1, klh.At(2, 3)*t11*t21+klv.At(2, 3)*t11*t21)
	kAdd(uy1, phiz1, klh.At(2, 3)*t11*t22+klv.At(2, 3)*t12*t21)

	kAdd(uz1, uz1, klh.At(2, 2)*t12*t12+klv.At(2, 2)*t22*t22)
	kAdd(uz1, phix1, klh.At(2, 3)*t12*t20+klv.At(2, 3)*t10*t22)
	kAdd(uz1, phiy1, klh.At(2, 3)*t12*t21+klv.At(2, 3)*t11*t22)
	kAdd(uz1, phiz1, klh.At(2, 3)*t12*t22+klv.At(2, 3)*t12*t22)

	kAdd(phix1, phix1, klh.At(3, 3)*t20*t20+klv.At(3, 3)*t10*t10)
	kAdd(phix1, phiy1, klh.At(3, 3)*t20*t21+klv.At(3, 3)*t10*t11)
	kAdd(phix1, phiz1, klh.At(3, 3)*t20*t22+klv.At(3, 3)*t10*t12)

	kAdd(phiy1, phiy1, klh.At(3, 3)*t21*t21+klv.At(3, 3)*t11*t11)
	kAdd(phiy1, phiz1, klh.At(3, 3)*t21*t22+klv.At(3, 3)*t11*t12)

	kAdd(phiz1, phiz1, klh.At(3, 3)*t22*t22+klv.At(3, 3)*t12*t12)

	rAdd(ux0, rlh.AtVec(0)*t10+rlv.AtVec(0)*t20)
	rAdd(uy0, rlh.AtVec(0)*t11+rlv.AtVec(0)*t21)
	rAdd(uz0, rlh.AtVec(0)*t12+rlv.AtVec(0)*t22)
	rAdd(phix0, rlh.AtVec(1)*t20+rlv.AtVec(1)*t10)
	rAdd(phiy0, rlh.AtVec(1)*t21+rlv.AtVec(1)*t11)
	rAdd(phiz0, rlh.AtVec(1)*t22+rlv.AtVec(1)*t12)
	rAdd(ux1, rlh.AtVec(2)*t10+rlv.AtVec(2)*t20)
	rAdd(uy1, rlh.AtVec(2)*t11+rlv.AtVec(2)*t21)
	rAdd(uz1, rlh.AtVec(2)*t12+rlv.AtVec(2)*t22)
	rAdd(phix1, rlh.AtVec(3)*t20+rlv.AtVec(3)*t10)
	rAdd(phiy1, rlh.AtVec(3)*t21+rlv.AtVec(3)*t11)
	rAdd(phiz1, rlh.AtVec(3)*t22+rlv.AtVec(3)*t12)
}

func (b *beam3d) indicesAsArray() *[12]Index {
	indices := [...]Index{
		{NodalID: b.n0.ID, Dof: Ux},
		{NodalID: b.n0.ID, Dof: Uy},
		{NodalID: b.n0.ID, Dof: Uz},
		{NodalID: b.n0.ID, Dof: Phix},
		{NodalID: b.n0.ID, Dof: Phiy},
		{NodalID: b.n0.ID, Dof: Phiz},
		{NodalID: b.n1.ID, Dof: Ux},
		{NodalID: b.n1.ID, Dof: Uy},
		{NodalID: b.n1.ID, Dof: Uz},
		{NodalID: b.n1.ID, Dof: Phix},
		{NodalID: b.n1.ID, Dof: Phiy},
		{NodalID: b.n1.ID, Dof: Phiz},
	}

	return &indices
}

func (b *beam3d) AddLoad(bc NeumannElementBC) bool {
	supported := false

	loadDispatch(bc,
		func(l *neumannConcentrated) {
			switch l.kind {
			case Uz, Uy, Phiy, Phiz:
				supported = true
			}
		},
		func(l *neumannConstant) {
			switch l.kind {
			case Uz, Uy:
				supported = true
			}
		},
		func(l *neumannLinear) {
			switch l.kind {
			case Uz, Uy:
				supported = true
			}
		},
	)

	return supported && b.oneDimElement.AddLoad(bc)
}

func (b *beam3d) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	return nil
}

func (b *beam3d) Indices(set map[Index]struct{}) {
	for _, index := range b.indicesAsArray() {
		set[index] = struct{}{}
	}
}
