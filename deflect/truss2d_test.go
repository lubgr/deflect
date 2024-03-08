package deflect

import (
	"maps"
	"math"
	"testing"

	"gonum.org/v1/gonum/mat"
)

var exampleMat = Material{
	CrossSection: &rectangular{0.1, 0.1},
	LinearElastic: LinearElastic{
		YoungsModulus: 100000.0 / (1e-3 * 1e-3),
		PoissonsRatio: 0.3,
		Density:       1,
	},
}

func truss2dTestNodes() (n0, n1 *Node) {
	n0 = &Node{ID: "A", X: 0, Z: 0}
	n1 = &Node{ID: "B", X: 1, Z: 0}
	return
}

func TestTruss2dCtorFailsOnZeroLength(t *testing.T) {
	n0 := &Node{ID: "A", X: 2, Y: 3, Z: -5}
	n1 := &Node{ID: "B", X: 2, Y: 3, Z: -5}

	actual, err := NewTruss2d("AB", n0, n1, &exampleMat, Truss2dDisjointNone)

	if err == nil || actual != nil {
		t.Errorf("Truss creation succeeded despite zero length, got %v", actual)
	}
}

func TestTruss2dCtorFailsOnEmptyElementId(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	actual, err := NewTruss2d("", n0, n1, &exampleMat, Truss2dDisjointNone)

	if err == nil || actual != nil {
		t.Errorf("Truss creation succeeded despite empty element id, got %v", actual)
	}
}

func TestTruss2dCtorSuccess(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	_, err := NewTruss2d("ID", n0, n1, &exampleMat, Truss2dDisjointNone)

	if err != nil {
		t.Errorf("Expected successful NewTruss2d call, got %v", err)
	}
}

func TestTruss2dReturnsID(t *testing.T) {
	// Stupid test, candidate for removal from day one...
	n0, n1 := truss2dTestNodes()
	id := "ID123"
	truss, _ := NewTruss2d(id, n0, n1, &exampleMat, Truss2dDisjointNone)

	if actual := truss.ID(); actual != id {
		t.Errorf("Expected truss to return an ID as %v, but got %v", id, actual)
	}
}

func TestTruss2dNumNodes(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	truss, _ := NewTruss2d("ID", n0, n1, &exampleMat, Truss2dDisjointNone)
	actual := truss.NumNodes()

	if actual != 2 {
		t.Errorf("Expected 2d truss to indicate 2 connected nodes, got %v", actual)
	}
}

func TestTruss2dIndexSet(t *testing.T) {
	n0, n1 := truss2dTestNodes()
	truss, _ := NewTruss2d("ID", n0, n1, &exampleMat, Truss2dDisjointNone)
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

func TestTruss2dRotationMatrix(t *testing.T) {
	// The rotation matrix is tested by applying the rotation to a point in the element's local
	// coordinate system and comparing the result to expected values.
	local := func(x, z float64) *mat.VecDense {
		return mat.NewVecDense(4, []float64{x, z, x, z})
	}
	global := local

	cases := [...]struct {
		x1, z1 float64
		local  *mat.VecDense
		global *mat.VecDense
	}{
		// The truss's first node is fixed at {2, 0, 2}, see below. All y coordinates are zero, given
		// the 2d domain.
		//
		// 1. Horizontal element
		{x1: 3, z1: 2, local: local(0, 0), global: global(0, 0)},
		{x1: 3, z1: 2, local: local(3, 2), global: global(3, -2)},
		{x1: 1, z1: 2, local: local(3, 2), global: global(-3, 2)},
		{x1: 1, z1: 2, local: local(3, -2), global: global(-3, -2)},
		// 2. Element rotated 45 degree counter clockwise
		{x1: 3, z1: 3, local: local(1, 0), global: global(1/math.Sqrt2, 1/math.Sqrt2)},
		{x1: 3, z1: 3, local: local(1, -1), global: global(0, math.Sqrt2)},
		{x1: 3, z1: 3, local: local(-1, -1), global: global(-math.Sqrt2, 0)},
		{x1: 3, z1: 3, local: local(1, 1), global: global(math.Sqrt2, 0)},
		// 3. Vertical element
		{x1: 2, z1: 10, local: local(1, 0), global: global(0, 1)},
		{x1: 2, z1: 10, local: local(1, 1), global: global(1, 1)},
		{x1: 2, z1: 10, local: local(-1, 1), global: global(1, -1)},
		{x1: 2, z1: -10, local: local(1, 1), global: global(-1, -1)},
		{x1: 2, z1: -10, local: local(5, 5), global: global(-5, -5)},
		{x1: 2, z1: -10, local: local(-5, 5), global: global(-5, 5)},
	}

	for _, test := range cases {
		n0 := &Node{ID: "A", X: 2, Z: 2}
		n1 := &Node{ID: "B", X: test.x1, Z: test.z1}
		truss := &truss2d{"ID", n0, n1, &exampleMat, Truss2dDisjointNone}
		actual := mat.NewVecDense(4, nil)
		actual.MulVec(truss.rotation(), test.local)

		if !mat.EqualApprox(test.global, actual, 1e-10) {
			t.Errorf("Expected point %v rotated to global XYZ to be %v, got %v",
				mat.Formatted(test.local.T()), mat.Formatted(test.global.T()), mat.Formatted(actual.T()))
		}
	}
}

func TestTruss2dLocalTangent(t *testing.T) {
	// Testing the local tangent is simpler than testing the global one. The latter is better captured
	// with a more integrated test that solves simple boundary value problems.
	tests := [...]struct {
		youngs, b, h  float64
		disjoint      Truss2dDisjoint
		k00, k02, k22 float64
	}{
		{youngs: 2, b: 1, h: 10, disjoint: Truss2dDisjointNone, k00: 4, k02: -4, k22: 4},
		{youngs: 3, b: 5, h: 1, disjoint: Truss2dDisjointNone, k00: 3, k02: -3, k22: 3},
		{youngs: 1, b: 2, h: 0.5, disjoint: Truss2dDisjointFirstNode},
		{youngs: 1, b: 1, h: 1, disjoint: Truss2dDisjointSecondNode},
	}

	for _, test := range tests {
		// Nodal coordinates are arranged so that the element length is 5
		n0 := &Node{ID: "A", X: 13, Z: -5}
		n1 := &Node{ID: "B", X: 10, Z: -1}
		material := exampleMat
		material.YoungsModulus = test.youngs
		material.CrossSection, _ = NewRectangularCrossSection(test.b, test.h)
		truss := &truss2d{"AB", n0, n1, &material, test.disjoint}
		actual := truss.localTangent()

		expected := mat.NewSymDense(4, nil)
		expected.SetSym(0, 0, test.k00)
		expected.SetSym(0, 2, test.k02)
		expected.SetSym(2, 2, test.k22)

		if !mat.EqualApprox(actual, expected, 1e-10) {
			t.Errorf("Expected local tangent to be \n%v\n, instead got \n%v\n",
				mat.Formatted(expected), mat.Formatted(actual))
		}
	}
}
