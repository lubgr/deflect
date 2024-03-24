package deflect

import (
	"fmt"

	"gonum.org/v1/gonum/mat"
)

// condenser instances apply static condensation on the element level, e.g., for hinged connections
// of truss or frame elements.
type condenser interface {
	// reduce eliminates degrees of freedoms and mutates the remaining entries in k and r to account
	// for the eliminated ones.
	reduce(k *mat.SymDense, r *mat.VecDense)
	// enhance uses k and r to restore the entry/entries in d that are disjoint through static
	// condensation.
	enhance(k *mat.SymDense, r, d *mat.VecDense)
}

// These types (except the noop noHinge) contain semi-generated implementations, and we don't expect
// a wild growth of more of these combinations, so isolating the cases through separate types is
// preferred over other dispatch mechanics. The main feature of these types if that they (1)
// pre-compute the static condensation for all needed scenarios, avoiding a matrix inversion and
// related failure handling at runtime, and (2) are not directly baked into the element
// formulations, so that material law and CO transformation from local to global axes is orthogonal
// to static condensation. One would typically assemble a vanilla (un-condensed) local element
// matrix, then apply static condensation through the condenser interface (which dispatches to the
// implementations below), and finally transform the matrices into global coordinates - static
// condensation is implemented so that it knows nothing of the other steps.
type noHinge struct{}
type preComputedHinge2x0 struct{}
type preComputedHinge2x1 struct{}
type preComputedHinge4x0 struct{}
type preComputedHinge4x1 struct{}
type preComputedHinge4x2 struct{}
type preComputedHinge4x3 struct{}
type preComputedHinge4x1x2 struct{}
type preComputedHinge4x1x3 struct{}
type preComputedHinge4x0x3 struct{}

func (*noHinge) reduce(k *mat.SymDense, r *mat.VecDense)     {}
func (*noHinge) enhance(k *mat.SymDense, r, d *mat.VecDense) {}

func (*preComputedHinge2x0) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k11 := k.At(0, 0), k.At(0, 1), k.At(1, 1)

	k.SetSym(0, 0, 0)
	k.SetSym(0, 1, 0)
	k.SetSym(1, 1, k11-k01*k01/k00)

	r0, r1 := r.AtVec(0), r.AtVec(1)

	r.SetVec(0, 0)
	r.SetVec(1, r1-k01*r0/k00)
}

func (*preComputedHinge2x0) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k00, k01 := k.At(0, 0), k.At(0, 1)
	r0 := r.AtVec(0)
	d1 := d.AtVec(1)

	d.SetVec(0, (r0-d1*k01)/k00)
}

func (*preComputedHinge2x1) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k11 := k.At(0, 0), k.At(0, 1), k.At(1, 1)

	k.SetSym(0, 0, k00-k01*k01/k11)
	k.SetSym(0, 1, 0)
	k.SetSym(1, 1, 0)

	r0, r1 := r.AtVec(0), r.AtVec(1)

	r.SetVec(0, r0-k01*r1/k11)
	r.SetVec(1, 0)
}

func (*preComputedHinge2x1) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k01, k11 := k.At(0, 1), k.At(1, 1)
	d0 := d.AtVec(0)
	r1 := r.AtVec(1)

	d.SetVec(1, (r1-d0*k01)/k11)
}

func (*preComputedHinge4x0) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(0, 0, 0)
	k.SetSym(0, 1, 0)
	k.SetSym(0, 2, 0)
	k.SetSym(0, 3, 0)
	k.SetSym(1, 1, k11-k01*k01/k00)
	k.SetSym(1, 2, k12-k01*k02/k00)
	k.SetSym(1, 3, k13-k01*k03/k00)
	k.SetSym(2, 2, k22-k02*k02/k00)
	k.SetSym(2, 3, k23-k02*k03/k00)
	k.SetSym(3, 3, k33-k03*k03/k00)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(0, 0)
	r.SetVec(1, r1-k01*r0/k00)
	r.SetVec(2, r2-k02*r0/k00)
	r.SetVec(3, r3-k03*r0/k00)

}
func (*preComputedHinge4x0) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	d1, d2, d3 := d.AtVec(1), d.AtVec(2), d.AtVec(3)
	r0 := r.AtVec(0)

	d.SetVec(0, (r0-d1*k01-d2*k02-d3*k03)/k00)
}

func (*preComputedHinge4x1) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(0, 0, k00-k01*k01/k11)
	k.SetSym(0, 1, 0)
	k.SetSym(0, 2, k02-k01*k12/k11)
	k.SetSym(0, 3, k03-k01*k13/k11)
	k.SetSym(1, 1, 0)
	k.SetSym(1, 2, 0)
	k.SetSym(1, 3, 0)
	k.SetSym(2, 2, k22-k12*k12/k11)
	k.SetSym(2, 3, k23-k12*k13/k11)
	k.SetSym(3, 3, k33-k13*k13/k11)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(0, r0-k01*r1/k11)
	r.SetVec(1, 0)
	r.SetVec(2, r2-k12*r1/k11)
	r.SetVec(3, r3-k13*r1/k11)
}

func (*preComputedHinge4x1) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k01 := k.At(0, 1)
	d0, d2, d3 := d.AtVec(0), d.AtVec(2), d.AtVec(3)
	r1 := r.AtVec(1)

	d.SetVec(1, (r1-d0*k01-d2*k12-d3*k13)/k11)
}

func (*preComputedHinge4x2) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(0, 0, k00-k02*k02/k22)
	k.SetSym(0, 1, k01-k02*k12/k22)
	k.SetSym(0, 2, 0)
	k.SetSym(0, 3, k03-k02*k23/k22)
	k.SetSym(1, 1, k11-k12*k12/k22)
	k.SetSym(1, 2, 0)
	k.SetSym(1, 3, k13-k12*k23/k22)
	k.SetSym(2, 2, 0)
	k.SetSym(2, 3, 0)
	k.SetSym(3, 3, k33-k23*k23/k22)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(0, r0-k02*r2/k22)
	r.SetVec(1, r1-k12*r2/k22)
	r.SetVec(2, 0)
	r.SetVec(3, r3-k23*r2/k22)
}

func (*preComputedHinge4x2) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k02, k12, k23, k22 := k.At(0, 2), k.At(1, 2), k.At(2, 3), k.At(2, 2)
	d0, d1, d3 := d.AtVec(0), d.AtVec(1), d.AtVec(3)
	r2 := r.AtVec(2)

	d.SetVec(2, (r2-d0*k02-d1*k12-d3*k23)/k22)
}

func (*preComputedHinge4x3) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(0, 0, k00-k03*k03/k33)
	k.SetSym(0, 1, k01-k03*k13/k33)
	k.SetSym(0, 2, k02-k03*k23/k33)
	k.SetSym(0, 3, 0)
	k.SetSym(1, 1, k11-k13*k13/k33)
	k.SetSym(1, 2, k12-k13*k23/k33)
	k.SetSym(1, 3, 0)
	k.SetSym(2, 2, k22-k23*k23/k33)
	k.SetSym(2, 3, 0)
	k.SetSym(3, 3, 0)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(0, r0-k03*r3/k33)
	r.SetVec(1, r1-k13*r3/k33)
	r.SetVec(2, r2-k23*r3/k33)
	r.SetVec(3, 0)
}

func (*preComputedHinge4x3) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k03, k13, k23, k33 := k.At(0, 3), k.At(1, 3), k.At(2, 3), k.At(3, 3)
	d0, d1, d2 := d.AtVec(0), d.AtVec(1), d.AtVec(2)
	r3 := r.AtVec(3)

	d.SetVec(3, (r3-d0*k03-d1*k13-d2*k23)/k33)
}

func (*preComputedHinge4x1x2) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(
		0,
		0,
		1.0/(k11*k22-(k12*k12))*(k00*(k11*k22-(k12*k12))-k01*(k01*k22-k02*k12)+k02*(k01*k12-k02*k11)),
	)
	k.SetSym(0, 1, 0)
	k.SetSym(0, 2, 0)
	k.SetSym(
		0,
		3,
		1.0/(k11*k22-(k12*k12))*(k03*(k11*k22-(k12*k12))-k13*(k01*k22-k02*k12)+k23*(k01*k12-k02*k11)),
	)
	k.SetSym(1, 1, 0)
	k.SetSym(1, 2, 0)
	k.SetSym(1, 3, 0)
	k.SetSym(2, 2, 0)
	k.SetSym(2, 3, 0)
	k.SetSym(
		3,
		3,
		1.0/(k11*k22-(k12*k12))*(k13*(k12*k23-k13*k22)-k23*(k11*k23-k12*k13)+k33*(k11*k22-(k12*k12))),
	)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(
		0,
		1.0/(k11*k22-(k12*k12))*(r0*(k11*k22-(k12*k12))-r1*(k01*k22-k02*k12)+r2*(k01*k12-k02*k11)),
	)
	r.SetVec(1, 0)
	r.SetVec(2, 0)
	r.SetVec(
		3,
		1.0/(k11*k22-(k12*k12))*(r1*(k12*k23-k13*k22)-r2*(k11*k23-k12*k13)+r3*(k11*k22-(k12*k12))),
	)
}

func (*preComputedHinge4x1x2) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k01, k02 := k.At(0, 1), k.At(0, 2)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	d0, d3 := d.AtVec(0), d.AtVec(3)
	r1, r2 := r.AtVec(1), r.AtVec(2)

	d.SetVec(1, (k12*(d0*k02+d3*k23-r2)-k22*(d0*k01+d3*k13-r1))*1.0/(k11*k22-(k12*k12)))
	d.SetVec(2, (-k11*(d0*k02+d3*k23-r2)+k12*(d0*k01+d3*k13-r1))*1.0/(k11*k22-(k12*k12)))
}

func (*preComputedHinge4x1x3) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(
		0,
		0,
		1.0/(k11*k33-(k13*k13))*(k00*(k11*k33-(k13*k13))-k01*(k01*k33-k03*k13)+k03*(k01*k13-k03*k11)),
	)
	k.SetSym(0, 1, 0)
	k.SetSym(
		0,
		2,
		1.0/(k11*k33-(k13*k13))*(k02*(k11*k33-(k13*k13))-k12*(k01*k33-k03*k13)+k23*(k01*k13-k03*k11)),
	)
	k.SetSym(0, 3, 0)
	k.SetSym(1, 1, 0)
	k.SetSym(1, 2, 0)
	k.SetSym(1, 3, 0)
	k.SetSym(
		2,
		2,
		1.0/(k11*k33-(k13*k13))*(-k12*(k12*k33-k13*k23)+k22*(k11*k33-(k13*k13))-k23*(k11*k23-k12*k13)),
	)
	k.SetSym(2, 3, 0)
	k.SetSym(3, 3, 0)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(
		0,
		1.0/(k11*k33-(k13*k13))*(r0*(k11*k33-(k13*k13))-r1*(k01*k33-k03*k13)+r3*(k01*k13-k03*k11)),
	)
	r.SetVec(1, 0)
	r.SetVec(
		2,
		1.0/(k11*k33-(k13*k13))*(-r1*(k12*k33-k13*k23)+r2*(k11*k33-(k13*k13))-r3*(k11*k23-k12*k13)),
	)
	r.SetVec(3, 0)

}

func (*preComputedHinge4x1x3) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k01, k03 := k.At(0, 1), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k23, k33 := k.At(2, 3), k.At(3, 3)
	d0, d2 := d.AtVec(0), d.AtVec(2)
	r1, r3 := r.AtVec(1), r.AtVec(3)

	d.SetVec(1, (k13*(d0*k03+d2*k23-r3)-k33*(d0*k01+d2*k12-r1))*1.0/(k11*k33-(k13*k13)))
	d.SetVec(3, (-k11*(d0*k03+d2*k23-r3)+k13*(d0*k01+d2*k12-r1))*1.0/(k11*k33-(k13*k13)))
}

func (*preComputedHinge4x0x3) reduce(k *mat.SymDense, r *mat.VecDense) {
	k00, k01, k02, k03 := k.At(0, 0), k.At(0, 1), k.At(0, 2), k.At(0, 3)
	k11, k12, k13 := k.At(1, 1), k.At(1, 2), k.At(1, 3)
	k22, k23 := k.At(2, 2), k.At(2, 3)
	k33 := k.At(3, 3)

	k.SetSym(0, 0, 0)
	k.SetSym(0, 1, 0)
	k.SetSym(0, 2, 0)
	k.SetSym(0, 3, 0)
	k.SetSym(
		1,
		1,
		1.0/(k00*k33-(k03*k03))*(-k01*(k01*k33-k03*k13)+k11*(k00*k33-(k03*k03))-k13*(k00*k13-k01*k03)),
	)
	k.SetSym(
		1,
		2,
		1.0/(k00*k33-(k03*k03))*(-k02*(k01*k33-k03*k13)+k12*(k00*k33-(k03*k03))-k23*(k00*k13-k01*k03)),
	)
	k.SetSym(1, 3, 0)
	k.SetSym(
		2,
		2,
		1.0/(k00*k33-(k03*k03))*(-k02*(k02*k33-k03*k23)+k22*(k00*k33-(k03*k03))-k23*(k00*k23-k02*k03)),
	)
	k.SetSym(2, 3, 0)
	k.SetSym(3, 3, 0)

	r0, r1, r2, r3 := r.AtVec(0), r.AtVec(1), r.AtVec(2), r.AtVec(3)

	r.SetVec(0, 0)
	r.SetVec(
		1,
		1.0/(k00*k33-(k03*k03))*(-r0*(k01*k33-k03*k13)+r1*(k00*k33-(k03*k03))-r3*(k00*k13-k01*k03)),
	)
	r.SetVec(
		2,
		1.0/(k00*k33-(k03*k03))*(-r0*(k02*k33-k03*k23)+r2*(k00*k33-(k03*k03))-r3*(k00*k23-k02*k03)),
	)
	r.SetVec(3, 0)

}

func (*preComputedHinge4x0x3) enhance(k *mat.SymDense, r, d *mat.VecDense) {
	k00, k02 := k.At(0, 0), k.At(0, 2)
	k01, k03 := k.At(0, 1), k.At(0, 3)
	k13, k23, k33 := k.At(1, 3), k.At(2, 3), k.At(3, 3)
	d1, d2 := d.AtVec(1), d.AtVec(2)
	r0, r3 := r.AtVec(0), r.AtVec(3)

	d.SetVec(0, (k03*(d1*k13+d2*k23-r3)-k33*(d1*k01+d2*k02-r0))*1.0/(k00*k33-(k03*k03)))
	d.SetVec(3, (-k00*(d1*k13+d2*k23-r3)+k03*(d1*k01+d2*k02-r0))*1.0/(k00*k33-(k03*k03)))
}

func newStaticHingeCondenser(all []Index, hinged map[Index]struct{}) (c condenser, err error) {
	indices := make([]int, 0, 2)

	for i, index := range all {
		if _, ok := hinged[index]; ok {
			indices = append(indices, i)
		}
	}

	if len(indices) == 0 {
		c = &noHinge{}
		return c, err
	}

	single := func(size, h int) bool {
		return len(all) == size && len(indices) == 1 && indices[0] == h
	}
	double := func(size, h0, h1 int) bool {
		return len(all) == size && len(indices) == 2 && indices[0] == h0 && indices[1] == h1
	}

	switch {
	case single(2, 0):
		c = &preComputedHinge2x0{}
	case single(2, 1):
		c = &preComputedHinge2x1{}
	case single(4, 0):
		c = &preComputedHinge4x0{}
	case single(4, 1):
		c = &preComputedHinge4x1{}
	case single(4, 2):
		c = &preComputedHinge4x2{}
	case single(4, 3):
		c = &preComputedHinge4x3{}
	case double(4, 1, 2):
		c = &preComputedHinge4x1x2{}
	case double(4, 1, 3):
		c = &preComputedHinge4x1x3{}
	case double(4, 0, 3):
		c = &preComputedHinge4x0x3{}
	default:
		err = fmt.Errorf("unsupported hinge setup, indices %v with hinges %v", all, indices)
	}

	return c, err
}
