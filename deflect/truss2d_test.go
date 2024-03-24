package deflect

import (
	"maps"
	"testing"

	"gonum.org/v1/gonum/mat"
)

func truss2dTestNodes() (n0, n1 *Node) {
	n0 = &Node{ID: "A", X: 0, Z: 0}
	n1 = &Node{ID: "B", X: 1, Z: 0}
	return
}

func TestTruss2dCtorSuccess(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	_, err := NewTruss2d("ID", n0, n1, &exampleMat, map[Index]struct{}{})

	if err != nil {
		t.Errorf("Expected successful NewTruss2d call, got %v", err)
	}
}

func TestTruss2dIndexSet(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	truss, _ := NewTruss2d("ID", n0, n1, &exampleMat, map[Index]struct{}{})
	actual := make(map[Index]struct{})
	expected := map[Index]struct{}{
		{NodalID: "A", Dof: Ux}: {},
		{NodalID: "A", Dof: Uz}: {},
		{NodalID: "B", Dof: Ux}: {},
		{NodalID: "B", Dof: Uz}: {},
	}

	truss.Indices(actual)

	if !maps.Equal(expected, actual) {
		t.Errorf("2d truss local dof: %v, expected %v", actual, expected)
	}
}

func TestTruss2dLocalTangent(t *testing.T) {
	// Testing the local tangent is simpler than testing the global one. The latter is better captured
	// with a more integrated test that solves simple boundary value problems.
	tests := [...]struct {
		youngs, b, h  float64
		k00, k01, k11 float64
	}{
		{youngs: 2, b: 1, h: 10, k00: 4, k01: -4, k11: 4},
		{youngs: 3, b: 5, h: 1, k00: 3, k01: -3, k11: 3},
	}

	for _, test := range tests {
		// Nodal coordinates are arranged so that the element length is 5
		n0 := &Node{ID: "A", X: 13, Z: -5}
		n1 := &Node{ID: "B", X: 10, Z: -1}
		material := exampleMat
		material.YoungsModulus = test.youngs
		material.CrossSection, _ = NewRectangularCrossSection(test.b, test.h)
		common, _ := newOneDimElement("AB", n0, n1, &material)
		truss := &truss2d{oneDimElement: common, hinges: nil}
		actual := truss.localNoHingeTangent(5)

		expected := mat.NewSymDense(2, nil)
		expected.SetSym(0, 0, test.k00)
		expected.SetSym(0, 1, test.k01)
		expected.SetSym(1, 1, test.k11)

		if !mat.EqualApprox(actual, expected, 1e-10) {
			t.Errorf("Expected local tangent to be \n%v\n, instead got \n%v\n",
				mat.Formatted(expected), mat.Formatted(actual))
		}
	}
}
