package elmt

import (
	"math"
	"testing"

	"github.com/lubgr/deflect/xyz"
)

func TestZeroLength(t *testing.T) {
	n0 := &xyz.Node{ID: "A", X: 2.0, Y: 3.0, Z: -5.0}
	n1 := &xyz.Node{ID: "B", X: 2.0, Y: 3.0, Z: -5.0}
	actual := length(n0, n1)

	if math.Abs(actual) > 0.e-20 {
		t.Errorf("Expected zero length, got %v", actual)
	}
}

func TestNonZeroLength(t *testing.T) {
	const expected = 24.351591323771842
	n0 := &xyz.Node{ID: "A", X: 2.0, Y: 10.0, Z: -5.0}
	n1 := &xyz.Node{ID: "B", X: -10.0, Y: 3.0, Z: 15.0}
	actual := length(n0, n1)

	if math.Abs(actual-expected) > 0.e-20 {
		t.Errorf("Expected length %v, got %v", expected, actual)
	}
}
