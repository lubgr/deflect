package bvp

import (
	"testing"

	"gonum.org/v1/gonum/mat"
)

func TestCholeskySolveSmallSuccessful(t *testing.T) {
	// Adjusted from gonum's cholesky_example_test.go
	A := mat.NewSymDense(4, []float64{
		120, 114, -4, -16,
		114, 118, 11, -24,
		-4, 11, 58, 17,
		-16, -24, 17, 73})

	b := mat.NewVecDense(4, []float64{1, 2, 3, 4})
	x := mat.NewVecDense(4, nil)

	ch := NewCholeskySolver()

	err := ch.SolveLinearSystem(A, b, x)

	if err != nil {
		t.Fatalf("Expected Solve to succeed, got %v", err)
	}

	expected := mat.NewVecDense(4, []float64{
		-0.239044214697,
		0.273229506199,
		-0.046809078910,
		0.103130891166})

	if !mat.EqualApprox(expected, x, 1e-10) {
		t.Errorf("Expected solution vector\n%v\nbut got\n%v", mat.Formatted(expected), mat.Formatted(x))
	}
}

func TestCholeskySolveSingularMatrix(t *testing.T) {
	dim := 150
	A := mat.NewSymDense(dim, nil)
	b := mat.NewVecDense(dim, nil)
	x := mat.NewVecDense(dim, nil)

	for i := range dim {
		ival := float64(i + 1)
		for j := i; j < dim; j++ {
			jval := float64(j + 1)
			A.SetSym(i, j, ival+jval)
		}
	}

	if mat.Det(A) != 0 {
		t.Fatalf("Input matrix was supposed to be singular")
	}

	ch := NewCholeskySolver()

	err := ch.SolveLinearSystem(A, b, x)

	if err == nil {
		t.Errorf("Expected Cholesky Ax=b solution to fail with singular A")
	}
}
