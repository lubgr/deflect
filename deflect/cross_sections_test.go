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

func TestRectangularIyy(t *testing.T) {
	b, h := 20.0, 30.0
	r, _ := NewRectangularCrossSection(b, h)
	expected := b * math.Pow(h, 3) / 12.0
	actual := r.Iyy()

	if math.Abs(actual-expected) > 1e-10 {
		t.Errorf("Expected Iyy to be %v, got %v", expected, actual)
	}
}

func TestRectangularNegativeDimensionsCauseConstructionError(t *testing.T) {
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

func TestConstantsCrossSectionConstruction(t *testing.T) {
	cases := []struct {
		params  map[string]float64
		failure bool
	}{
		{params: map[string]float64{"A": 1.0, "Iyy": 1.0}, failure: false},
		{params: map[string]float64{"A": -1.0, "Iyy": 1.0}, failure: true},
		{params: map[string]float64{"A": 1.0, "Iyy": -1.0}, failure: true},
		{params: map[string]float64{"A": 0, "Iyy": 0}, failure: true},
	}

	for _, test := range cases {
		_, err := NewConstantsCrossSections(test.params)

		if test.failure == (err == nil) {
			t.Errorf(
				"Unexpected constants creation result with %v, expected failure: %v",
				test.params,
				test.failure,
			)
		}
	}
}
