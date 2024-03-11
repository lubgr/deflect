package deflect

import (
	"math"
	"testing"
)

func TestZeroLength(t *testing.T) {
	n0 := &Node{ID: "A", X: 2.0, Y: 3.0, Z: -5.0}
	n1 := &Node{ID: "B", X: 2.0, Y: 3.0, Z: -5.0}
	actual := length(n0, n1)

	if math.Abs(actual) > 0.e-20 {
		t.Errorf("Expected zero length, got %v", actual)
	}
}

func TestNonZeroLength(t *testing.T) {
	const expected = 24.351591323771842
	n0 := &Node{ID: "A", X: 2.0, Y: 10.0, Z: -5.0}
	n1 := &Node{ID: "B", X: -10.0, Y: 3.0, Z: 15.0}
	actual := length(n0, n1)

	if math.Abs(actual-expected) > 0.e-20 {
		t.Errorf("Expected length %v, got %v", expected, actual)
	}
}

func TestSineCosine2d(t *testing.T) {
	cases := []struct{ x1, z1, angle float64 }{
		{x1: 0, z1: 0, angle: 0},
		{x1: 1, z1: 0, angle: 0},
		{x1: 1, z1: 1, angle: 45},
		{x1: 0, z1: 1, angle: 90},
		{x1: -1, z1: 1, angle: 135},
		{x1: -1, z1: 0, angle: 180},
		{x1: -1, z1: -1, angle: 225},
		{x1: 0, z1: -1, angle: 270},
		{x1: math.Cos(1.23456), z1: math.Sin(1.23456), angle: 1.23456 * 180.0 / math.Pi},
	}

	for _, test := range cases {
		n0 := &Node{ID: "A", X: 0.0, Y: 0.0, Z: 0.0}
		n1 := &Node{ID: "B", X: test.x1, Y: 0.0, Z: test.z1}
		sin, cos := sineCosine2d(n0, n1)
		alpha := math.Pi * test.angle / 180.0

		if expected := math.Sin(alpha); math.Abs(sin-expected) > 1.e-10 {
			t.Errorf("Expected angle %v with sine %v, got sine %v", test.angle, expected, sin)
		}

		if expected := math.Cos(alpha); math.Abs(cos-expected) > 1.e-10 {
			t.Errorf("Expected angle %v with cosine %v, got cosine %v", test.angle, expected, cos)
		}
	}
}
