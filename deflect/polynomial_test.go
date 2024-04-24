package deflect

import (
	"cmp"
	"slices"
	"testing"

	"gonum.org/v1/gonum/floats/scalar"
)

func TestPolyPieceDegree(t *testing.T) {
	cases := []struct {
		coeff  []float64
		expect int
	}{
		{coeff: []float64{0}, expect: 0},
		{coeff: []float64{0, 1}, expect: 1},
		{coeff: []float64{10, 10, 10}, expect: 2},
		{coeff: []float64{10, 0, 0, 0, 0, 0, 0, -5}, expect: 7},
	}

	for _, test := range cases {
		p := PolyPiece{X0: 0, XE: 1, Coeff: test.coeff}
		actual := p.Degree()

		if actual != test.expect {
			t.Errorf("Expected degree of %v to be %v, got %v", p, test.expect, actual)
		}
	}
}

func TestPolySequenceIntervalIndex(t *testing.T) {
	haystack := PolySequence{
		{X0: 0, XE: 1, Coeff: nil},
		{X0: 1, XE: 5, Coeff: nil},
		{X0: 5, XE: 5.5, Coeff: nil},
		{X0: 5.5, XE: 7, Coeff: nil},
		{X0: 7, XE: 10, Coeff: nil},
	}
	cases := []struct {
		x0, xE float64
		expect int
	}{
		{x0: 0, xE: 1, expect: 0},
		{x0: 0.1, xE: 1, expect: -1},
		{x0: 0, xE: 0.99, expect: -1},
		{x0: 1, xE: 5, expect: 1},
		{x0: 5, xE: 5.5, expect: 2},
		{x0: 5, xE: 7, expect: -1},
		{x0: 7, xE: 10, expect: 4},
		{x0: 10, xE: 10, expect: -1},
	}

	for _, test := range cases {
		actual := haystack.IntervalIndex(test.x0, test.xE)

		if actual != test.expect {
			t.Errorf(
				"Expected %vth polynomial to be found for query interval [%v, %v], got %v",
				test.expect,
				test.x0,
				test.xE,
				actual,
			)
		}
	}
}

func TestPolyPieceEval(t *testing.T) {
	cases := []struct {
		x0, xE, x, expect float64
		coeff             []float64
		success           bool
	}{
		{x0: 0, xE: 1, coeff: []float64{0}, x: 0, expect: 0, success: true},
		{x0: -3, xE: 5, coeff: []float64{0}, x: -3.5, expect: 0, success: false},
		{x0: -3, xE: 5, coeff: []float64{0}, x: 10, expect: 0, success: false},
		{x0: 0, xE: 1, coeff: []float64{2, 3, 4, 5}, x: 0.5, expect: 5.125, success: true},
		{
			x0:      -2,
			xE:      4,
			coeff:   []float64{0, -1, 2.25, 0, 7, 0, -5},
			x:       -1.5,
			expect:  -14.953125,
			success: true,
		},
	}

	for _, test := range cases {
		p := PolyPiece{X0: test.x0, XE: test.xE, Coeff: test.coeff}
		actual, err := p.Eval(test.x)

		if test.success && err != nil {
			t.Errorf("Unexpected failure: %v", err)
		} else if !scalar.EqualWithinAbsOrRel(actual, test.expect, 1e-8, 1e-8) {
			t.Errorf("Expected %v at x = %v to evaluate to %v, got %v", p, test.x, test.expect, actual)
		}
	}
}

func TestPolySequenceFlatten(t *testing.T) {
	cases := []struct {
		name          string
		input, expect PolySequence
	}{
		{
			name:   "Nil",
			input:  nil,
			expect: PolySequence{},
		},
		{
			name:   "Empty",
			input:  PolySequence{},
			expect: PolySequence{},
		},
		{
			name:   "Single",
			input:  PolySequence{{X0: 0, XE: 5, Coeff: []float64{0, 0, 3}}},
			expect: PolySequence{{X0: 0, XE: 5, Coeff: []float64{0, 0, 3}}},
		},
		{
			name: "Keep trailing zeros",
			input: PolySequence{
				{X0: -2, XE: 2, Coeff: []float64{1, 2, 3}},
				{X0: -2, XE: 2, Coeff: []float64{-1, -2, -3}},
			},
			expect: PolySequence{{X0: -2, XE: 2, Coeff: []float64{0, 0, 0}}},
		},
		{
			name: "Same domain",
			input: PolySequence{
				{X0: 2, XE: 4, Coeff: []float64{1, 2, 3}},
				{X0: 2, XE: 4, Coeff: []float64{2, 3, 4}},
				{X0: 2, XE: 4, Coeff: []float64{-3, 4, -1}},
				{X0: 2, XE: 4, Coeff: []float64{4, 0, 7}},
			},
			expect: PolySequence{
				{X0: 2, XE: 4, Coeff: []float64{4, 9, 13}},
			},
		},
		{
			name: "Overlapping 1",
			input: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1.5, Coeff: []float64{1, 2, 3}},
			},
			expect: PolySequence{
				{X0: 0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1, Coeff: []float64{2, 4, 6}},
				{X0: 1, XE: 1.5, Coeff: []float64{1, 2, 3}},
			},
		},
		{
			name: "Overlapping 2",
			input: PolySequence{
				{X0: 0, XE: 2, Coeff: []float64{1, 2, 3}},
				{X0: 1.25, XE: 1.9, Coeff: []float64{-1, -2, -4}},
				{X0: 1, XE: 1.25, Coeff: []float64{4}},
				{X0: 0, XE: 0.5, Coeff: []float64{0, 0, 3}},
				{X0: 0.5, XE: 1, Coeff: []float64{1, 0, 3}},
				{X0: 1, XE: 1.25, Coeff: []float64{5}},
				{X0: 0, XE: 1, Coeff: []float64{0, 1, -13}},
			},
			expect: PolySequence{
				{X0: 0, XE: 0.5, Coeff: []float64{1, 3, -7}},
				{X0: 0.5, XE: 1, Coeff: []float64{2, 3, -7}},
				{X0: 1, XE: 1.25, Coeff: []float64{10, 2, 3}},
				{X0: 1.25, XE: 1.9, Coeff: []float64{0, 0, -1}},
				{X0: 1.9, XE: 2, Coeff: []float64{1, 2, 3}},
			},
		},
		{
			name: "With gap",
			input: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1, Coeff: []float64{0, -2, 1}},
				{X0: 1.5, XE: 2, Coeff: []float64{4, 5, 6}},
			},
			expect: PolySequence{
				{X0: 0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1, Coeff: []float64{1, 0, 4}},
				{X0: 1.5, XE: 2, Coeff: []float64{4, 5, 6}},
			},
		},
	}

	for _, test := range cases {
		t.Run(test.name, func(t *testing.T) {
			actual := test.input.flatten()

			if !polyEqual(actual, test.expect, 1e-10) {
				t.Errorf("Expected polynomial %v to flatten to %v, got %v", test.input, test.expect, actual)
			}
		})
	}
}

func TestPolySequenceCompactIdentical(t *testing.T) {
	cases := []struct {
		ps, expect PolySequence // If expect is empty, the output must be identical to the ps
		tol        float64      // If zero-initialised, a default tolerance is used
	}{
		{
			ps:     PolySequence{},
			expect: PolySequence{},
		},
		{
			ps: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}},
			},
		},
		{
			ps: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}},
				{X0: 1, XE: 2, Coeff: []float64{1, 2}},
			},
		},
		{
			ps: PolySequence{
				{X0: 0.0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 2.0, Coeff: []float64{1, 2, 3}},
			},
			expect: PolySequence{
				{X0: 0, XE: 2, Coeff: []float64{1, 2, 3}},
			},
		},
		{
			ps: PolySequence{
				{X0: 0.0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 2.0, Coeff: []float64{1, 2 + 2e-10, 3}},
			},
		},
		{
			ps: PolySequence{
				{X0: 0.0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 2.0, Coeff: []float64{1, 2 + 2e-10, 3}},
			},
			expect: PolySequence{
				{X0: 0.0, XE: 2.0, Coeff: []float64{1, 2, 3}},
			},
			tol: 1e-9,
		},
		{
			ps: PolySequence{
				{X0: 0.0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1.0, Coeff: []float64{1, 2, 3}},
				{X0: 1.0, XE: 1.5, Coeff: []float64{1, 2, 3}},
				{X0: 1.5, XE: 2.0, Coeff: []float64{1, 2, 3}},
			},
			expect: PolySequence{
				{X0: 0, XE: 2, Coeff: []float64{1, 2, 3}},
			},
		},
		{
			ps: PolySequence{
				{X0: 0.0, XE: 0.5, Coeff: []float64{1, 2, 3}},
				{X0: 0.5, XE: 1.0, Coeff: []float64{1, 2, 3}},
				{X0: 1.0, XE: 1.5, Coeff: []float64{1, 2, 3, 4}},
				{X0: 1.5, XE: 1.5, Coeff: []float64{1, 2, 3}},
				{X0: 1.5, XE: 2.0, Coeff: []float64{1, 2, 3}},
				// With a gap in between, we don't join:
				{X0: 3.0, XE: 4.0, Coeff: []float64{0}},
				{X0: 4.0, XE: 5.0, Coeff: []float64{0}},
			},
			expect: PolySequence{
				{X0: 0.0, XE: 1.0, Coeff: []float64{1, 2, 3}},
				{X0: 1.0, XE: 1.5, Coeff: []float64{1, 2, 3, 4}},
				{X0: 1.5, XE: 2.0, Coeff: []float64{1, 2, 3}},
				{X0: 3.0, XE: 5.0, Coeff: []float64{0}},
			},
		},
	}

	for _, test := range cases {
		tol := cmp.Or(test.tol, 1e-10)
		orig := slices.Clone(test.ps)
		test.ps = test.ps.CompactIdentical(tol)

		if len(test.expect) == 0 {
			test.expect = orig
		}

		if !polyEqual(test.ps, test.expect, 1e-10) {
			t.Errorf("Expected %v to be compacted to %v, got %v", orig, test.expect, test.ps)
		}
	}
}

func TestPolySequenceTrimTrailingZeros(t *testing.T) {
	cases := []struct {
		ps, expect PolySequence // If expect is empty, the output must be identical to the input
		tol        float64      // If zero-initialised, a default tolerance is used
	}{
		{ps: PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}}}},
		{ps: PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3, 1e-8}}}},
		{
			ps:     PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3, 1e-8, 1e-8, 1e-8}}},
			expect: PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}}},
			tol:    5e-7,
		},
		{
			ps:     PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3, 0}}},
			expect: PolySequence{{X0: 0, XE: 1, Coeff: []float64{1, 2, 3}}},
		},
		{
			ps: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3, 0, 0, 4}},
				{X0: 1, XE: 2, Coeff: []float64{0, 0, 0, 0, 0}},
				{X0: 2, XE: 3, Coeff: []float64{1, 2, 3, 0, 0, -1e-10}},
			},
			expect: PolySequence{
				{X0: 0, XE: 1, Coeff: []float64{1, 2, 3, 0, 0, 4}},
				{X0: 1, XE: 2, Coeff: []float64{0}},
				{X0: 2, XE: 3, Coeff: []float64{1, 2, 3}},
			},
		},
	}

	for _, test := range cases {
		tol := cmp.Or(test.tol, 1e-10)
		orig := slices.Clone(test.ps)

		if len(test.expect) == 0 {
			test.expect = append(test.expect, test.ps...)
		}

		test.ps.TrimTrailingZeros(tol)

		if !polyEqual(test.ps, test.expect, 1e-10) {
			t.Errorf("Expected %v to be zero-trimmed to %v, got %v", orig, test.expect, test.ps)
		}
	}
}

// polyEqual returns true if p0 and p1 are approximately the same, using absolute and relative
// floating point comparison with the given tolerance.
func polyEqual(p0, p1 PolySequence, tol float64) bool {
	return slices.EqualFunc(p0, p1, func(a, b PolyPiece) bool { return polyDataEqual(a, b, tol) })
}

func polyDataEqual(p0, p1 PolyPiece, tol float64) bool {
	approxEq := func(a, b float64) bool {
		return scalar.EqualWithinAbsOrRel(a, b, tol, tol)
	}

	sameX0 := approxEq(p0.X0, p1.X0)
	sameXE := approxEq(p0.XE, p1.XE)
	sameCoeff := slices.EqualFunc(p0.Coeff, p1.Coeff, approxEq)

	return sameX0 && sameXE && sameCoeff
}
