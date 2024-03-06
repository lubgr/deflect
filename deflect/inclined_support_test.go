package deflect

import (
	"math"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestInclinedSupportTransformation(t *testing.T) {
	// This test makes sure the coordinate transformation implementation of an inclined support gives
	// identical results to computing the transformation with a full dense transformation matrix. The
	// latter is known to work correctly (using Gonum operations), while the former implements a
	// custom memory-efficient scheme and might change over time, hence the test to prevent any
	// regressions from these potential changes (it's arguably not a great test, too bulky and
	// dependent on the implementation details of the function under test).
	cases := []float64{0, 45, 90, 123.456}
	from, to := Index{NodalID: "A", Dof: Ux}, Index{NodalID: "A", Dof: Uy}
	indices := map[Index]int{from: 3, to: 20}
	dim := 30

	for _, angle := range cases {
		alpha := angle * math.Pi / 180.0
		transformer, _ := NewInclinedSupport(from, to, alpha)

		k, rhs, d := matricesToTransform(dim)
		kref, rhsref, dref := referenceMatricesToTransform(k, rhs, d)
		r := transformationMatrix(dim, 3, 20, alpha)
		kref.Mul(r.T(), kref)
		kref.Mul(kref, r)
		dref.MulVec(r.T(), d)
		rhsref.MulVec(r.T(), rhs)

		transformer.Pre(indices, k, rhs, d)

		if !mat.EqualApprox(k, kref, 1e-8) {
			t.Errorf("Expected Pre operation to compute tᵀ·k·t, but reference result differs")
		}
		if !mat.EqualApprox(d, dref, 1e-8) {
			t.Errorf("Expected Pre operation to compute tᵀ·d, but reference result differs")
		}

		transformer.Post(indices, rhs, d)

		dref.MulVec(r, dref)
		rhsref.MulVec(r, rhsref)

		if !mat.EqualApprox(d, dref, 1e-8) {
			t.Errorf("Expected Post to compute t·d, but reference result differs")
		}
		if !mat.EqualApprox(rhs, rhsref, 1e-8) {
			t.Errorf("Expected Post to compute t·rhs, but reference result differs")
		}
	}
}

func matricesToTransform(dim int) (k *mat.SymDense, rhs, d *mat.VecDense) {
	rhs, d = mat.NewVecDense(dim, nil), mat.NewVecDense(dim, nil)
	k = mat.NewSymDense(dim, nil)

	value := 0.0
	for i := range dim {
		rhs.SetVec(i, value)
		d.SetVec(i, 2*value)
		for j := i; j < dim; j++ {
			k.SetSym(i, j, value)
			value += 1.0
		}
	}

	return k, rhs, d
}

func referenceMatricesToTransform(
	k *mat.SymDense,
	rhs, d *mat.VecDense,
) (kref *mat.Dense, rhsref, dref *mat.VecDense) {
	dim := k.SymmetricDim()
	rhsref, dref = mat.NewVecDense(dim, nil), mat.NewVecDense(dim, nil)
	kref = mat.NewDense(dim, dim, nil)

	kref.Copy(k)
	dref.CopyVec(d)
	rhsref.CopyVec(rhs)

	return kref, rhsref, dref
}

func transformationMatrix(dim, k, l int, alpha float64) *mat.Dense {
	t := mat.NewDense(dim, dim, nil)

	for i := 0; i < dim; i++ {
		t.Set(i, i, 1.0)
	}

	t.Set(k, k, math.Cos(alpha))
	t.Set(k, l, -math.Sin(alpha))
	t.Set(l, l, math.Cos(alpha))
	t.Set(l, k, math.Sin(alpha))

	return t
}
