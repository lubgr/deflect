package deflect

import (
	"fmt"
	"slices"

	"gonum.org/v1/gonum/floats/scalar"
)

// PolyPiece describes a one-dimensional piecewise polynomial with non-negative exponents.
type PolyPiece struct {
	// Start and end of the domain. Precisely, the interval is [X0, XE], so that both values are
	// included in the domain.
	X0, XE float64
	// The index of each coefficient is the associated exponent. Coefficients can be zero, but never
	// the last one except when there is only a single entry. For example, 3·x² - 2.5 would be
	// represented as []float{-2.5, 0.0, 3.0}.
	Coeff []float64
}

// Degree returns the degree of the given polynomial.
func (p *PolyPiece) Degree() int {
	return len(p.Coeff) - 1
}

// Eval evaluates the polynomial at x. If x is outside of the domain, an error is returned.
func (p *PolyPiece) Eval(x float64) (float64, error) {
	if x+1e-10 < p.X0 || x-1e10 > p.XE {
		return 0, fmt.Errorf("%v outside of domain [%v, %v]", x, p.X0, p.XE)
	}

	y := p.Coeff[0]

	for xp, i := x, 1; i < len(p.Coeff); i++ {
		y += p.Coeff[i] * xp
		xp *= x
	}

	return y, nil
}

// PolySequence is a piecewise polynomial.
type PolySequence []PolyPiece

// IntervalIndex behaves like [slices.Index]. A match is determined by comparing x0 and xE with the
// domain of the polynomials in ps, where both x0 and xE must match.
func (ps PolySequence) IntervalIndex(x0, xE float64) int {
	absTol, relTol := 1e-10, 1e-10

	for i, p := range ps {
		equalX0 := scalar.EqualWithinAbsOrRel(x0, p.X0, absTol, relTol)
		equalXE := scalar.EqualWithinAbsOrRel(xE, p.XE, absTol, relTol)

		if equalX0 && equalXE {
			return i
		}
	}

	return -1
}

// flatten combines the given, additive piecewise polynomials so that there is only a single
// PolyPiece instance per sub-interval of the domain, and the sub-intervals do not overlap. The
// number of returned piecewise polynomials can be larger or smaller than the number of given
// polynomials, or identical. Example: given three polynomials over the domains [0, a], [0, a/3],
// and [a/3, 2a/3], the polynomials [0, a/3], [a/3, 2a/3], [2a/3, a] are returned, with summed
// polynomial coefficients where appropriate (e.g., the coefficients of the input polynomial over
// [0, a] are added to the output coefficients of all piecewise polynomials). It is possible that
// the resulting polynomial has trailing zero or near-zero coefficients, flatten does not try to
// trim these; in this sense, polynomials may be returned that are not canonical. We accept this
// since the approximate nature of the trimming step with an ideally numerical zero tolerance is
// inherently context-dependent. Use [PolySequence.TrimTrailingZeros] for this as a post-processing
// step, usually followed by [PolySequence.CompactIdentical].
func (ps PolySequence) flatten() PolySequence {
	xs := make([]float64, 0, 2*len(ps))
	result := make(PolySequence, 0, len(ps))

	for _, p := range ps {
		xs = append(xs, p.X0, p.XE)
	}

	slices.Sort(xs)

	approxEq := func(a, b float64) bool { return scalar.EqualWithinAbs(a, b, 1e-10) }
	sameRange := func(p PolyPiece, x0, xE float64) bool {
		return (p.X0 < x0 || approxEq(p.X0, x0)) && (p.XE > xE || approxEq(p.XE, xE))
	}

	xs = slices.CompactFunc(xs, approxEq)

	// If ps is not empty, xs must have at least 2 elements
	for i := 0; i < len(xs)-1; i++ {
		x0, xE := xs[i], xs[i+1]
		interval := PolyPiece{X0: x0, XE: xE, Coeff: nil}

		for _, p := range ps {
			if !sameRange(p, x0, xE) {
				continue
			}

			// Make space for new coefficients
			current := len(interval.Coeff)
			for n := current; n < max(current, len(p.Coeff)); n++ {
				interval.Coeff = append(interval.Coeff, 0.0)
			}

			for j, coeff := range p.Coeff {
				interval.Coeff[j] += coeff
			}
		}

		// Coeff can be nil with gaps between intervals. Should not happen when used for element
		// interpolations, but we should still handle it well.
		if interval.Coeff != nil {
			result = append(result, interval)
		}
	}

	return result
}

// CompactIdentical squashes consecutive identical polynomials into a single one, i.e. it mutates
// the object in place if there are consecutive identical PolyPiece instances. Since the size of the
// underlying slice may change, callers must use the return value the same way they would when
// appending to a slice.
func (ps PolySequence) CompactIdentical(tol float64) PolySequence {
	approxEq := func(a, b float64) bool { return scalar.EqualWithinAbsOrRel(a, b, tol, tol) }
	// Store the indices to polynomials that shall be removed. It would be easier to aggregate a new
	// slice with polynomials to keep instead, but the dominant scenario here is no compacting at all.
	var removal []int

	for i := 1; i < len(ps); i++ {
		if !approxEq(ps[i-1].XE, ps[i].X0) {
			// There is a gap between two subsequent polynomials, should not happen for element
			// interpolations, but let's treat this gracefully.
			continue
		}

		if slices.EqualFunc(ps[i-1].Coeff, ps[i].Coeff, approxEq) {
			removal = append(removal, i)
		}
	}

	// Probably no need to optimise this operation, removal will be empty > 99% of the time.
	for i := len(removal) - 1; i >= 0; i-- {
		r := removal[i]
		ps[r-1].XE = ps[r].XE
		ps = slices.Delete(ps, r, r+1)
	}

	return ps
}

// TrimTrailingZeros strips trailing approximate zeros in place, except the very first one, which is
// truncated if approximately zero.
func (ps PolySequence) TrimTrailingZeros(zeroTol float64) {
	approxZero := func(x float64) bool {
		return scalar.EqualWithinAbs(x, 0.0, zeroTol)
	}

	for i := 0; i < len(ps); i++ {
		for n := len(ps[i].Coeff) - 1; n > 0; n-- {
			if approxZero(ps[i].Coeff[n]) {
				ps[i].Coeff = ps[i].Coeff[:n]
			} else {
				break
			}
		}

		if approxZero(ps[i].Coeff[0]) {
			ps[i].Coeff[0] = 0
		}
	}
}
