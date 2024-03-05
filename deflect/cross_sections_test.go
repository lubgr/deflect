package deflect

import (
	"math"
	"testing"
)

func TestRectangularArea(t *testing.T) {
	r, _ := NewRectangularCrossSection(12.3, 45.6)
	expected := 12.3 * 45.6
	actual := r.Area()

	if math.Abs(actual-expected) > 1e-10 {
		t.Errorf("Expected area of rectangle to be %v, got %v", expected, actual)
	}
}

func TestNegativeDimensionsCauseConstructionError(t *testing.T) {
	cases := []struct{ b, h float64 }{
		{b: -1, h: 1},
		{b: 1, h: -1},
		{b: -1, h: -1},
	}

	for _, test := range cases {
		if _, err := NewRectangularCrossSection(test.b, test.h); err == nil {
			t.Errorf("Expected rectangular profile construction to fail with b/h = %v/%v", test.b, test.h)
		}
	}
}
