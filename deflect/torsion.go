package deflect

import (
	"errors"
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// NewTorsionBar returns a new linear-elastic 1d torsion element.
func NewTorsionBar(
	id string,
	n0, n1 *Node,
	material *Material,
	hinges map[Index]struct{},
) (Element, error) {
	localIndices := [...]Index{{NodalID: n0.ID, Dof: Phix}, {NodalID: n1.ID, Dof: Phix}}

	common, errCommon := newOneDimElement(id, n0, n1, material)
	condenser, errCondensation := newStaticHingeCondenser(localIndices[:], hinges)

	if err := errors.Join(errCommon, errCondensation); err != nil {
		return nil, fmt.Errorf("failed to instantiate new torsion bar: %w", err)
	}

	return &torsion{oneDimElement: common, hinges: condenser}, nil
}

type torsion struct {
	oneDimElement
	hinges condenser
}

func (t *torsion) Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense) {
	l := length(t.n0, t.n1)
	kl := t.localNoHingeTangent(l)
	fmt.Println(kl)
}

func (t *torsion) localNoHingeTangent(l float64) *mat.SymDense {
	k := mat.NewSymDense(2, nil)
	GIp := t.material.ShearModulus() * t.material.Ixx()

	k.SetSym(0, 0, GIp/l)
	k.SetSym(0, 1, -GIp/l)
	k.SetSym(1, 1, GIp/l)

	return k
}

func (t *torsion) indicesAsArray() *[2]Index {
	indices := [...]Index{
		{NodalID: t.n0.ID, Dof: Phix},
		{NodalID: t.n1.ID, Dof: Phix},
	}

	return &indices
}

func (t *torsion) AddLoad(bc NeumannElementBC) bool {
	var kind Dof

	loadDispatch(bc,
		func(load *neumannConcentrated) { kind = load.kind },
		func(load *neumannConstant) { kind = load.kind },
		func(load *neumannLinear) { kind = load.kind },
	)

	if kind == Phix {
		return t.oneDimElement.AddLoad(bc)
	}

	return false
}

func (t *torsion) Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence {
	switch which {
	case FctPhix, FctMx:
	default:
		return nil
	}

	return nil
}

func (t *torsion) Indices(set map[Index]struct{}) {
	for _, index := range t.indicesAsArray() {
		set[index] = struct{}{}
	}
}
