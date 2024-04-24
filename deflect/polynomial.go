package deflect

import (
	"fmt"
	"log"
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

// derive computes the polynomial's derivative and returns it as a new, independent object. Any
// constant polynomial (single coefficient in the zero slot, zero or non-zero) will lead to a
// zero-polynomial. The polynomial domain is copied as is.
func (p *PolyPiece) derive() PolyPiece {
	diff := PolyPiece{X0: p.X0, XE: p.XE, Coeff: nil}

	switch len(p.Coeff) {
	case 0:
		// Shouldn't happen in practice, but no harm to handle it gracefully with a no-op.
	case 1:
		diff.Coeff = []float64{0}
	default:
		// For example, with p = 3·x⁵ - 2·x² + x + 5, we have the coefficients
		//   p.Coeff = {5, 1, -2, 0, 0, 3}
		//   indices =  0  1   2  3  4  5
		// The derivative is computed by dropping the first coefficient, shifting everything "to the
		// left", and multiplying with the exponent:
		//   d/dx = {1·1, -2·2, 0·3, 0·4, 3·5} = {1, -4, 0, 0, 15}
		//   indices =													 {0,  1, 2, 3, 4}
		// which results in 15·x⁴ - 4·x + 1 by
		diff.Coeff = make([]float64, len(p.Coeff)-1)
		for j := 1; j < len(p.Coeff); j++ {
			diff.Coeff[j-1] = float64(j) * p.Coeff[j]
		}
	}

	return diff
}

// integrate does what the name suggests. An integration constant can be passed, which is used such
// that [PolyPiece.PolyEval] of p at p.X0 is constant. The polynomial domain is left as is. Note
// that integrating a near-zero constant polynomial, e.g. coefficients are {1e-8}, will result in a
// linear polynomial. integrate does not attempt to avoid this, but it can make sense to pass the
// result through [PolySequence.TrimTrailingZeros] after integration to drop the near-zero linear
// term.
func (p *PolyPiece) integrate(constant float64) PolyPiece {
	integrated := PolyPiece{X0: p.X0, XE: p.XE, Coeff: nil}

	// For example, with p = 3·x⁵ - 2·x² + x + 5, we have the coefficients
	//   p.Coeff = {5, 1, -2, 0, 0, 3}
	// The integral is computed by shifting everything "to the right", dividing by the incremented
	// exponent, and inserting the integration constant minus the evaluation of the integrated
	// function at X0:
	//   ∫ {5, 1, -2, 0, 0, 3} dx = {constant - ∫...dx | X0, 5/2, 1/3, -2/4, 0, 0, 3/6}.
	integrated.Coeff = make([]float64, len(p.Coeff)+1)

	for j := 0; j < len(p.Coeff); j++ {
		integrated.Coeff[j+1] = p.Coeff[j] / float64(j+1)
	}

	lower, errEval := integrated.Eval(integrated.X0)

	if errEval != nil {
		log.Printf("Bug: evaluating polynomial at lower boundary must not fail: %v", errEval)
	}

	integrated.Coeff[0] = constant - lower

	return integrated
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

// integrate uses [PolyPiece.integrate] to integrate individual polynomials. The receiver is
// unchanged, and the integral is returned as a new object. The sequence of polynomials is expected
// to have contiguous domains, i.e., gaps and overlaps must not exist between intervals. The
// integration constant is used for first polynomial in the sequence, all other integration
// constants are determined by evaluating the precedent polynomial.
func (ps PolySequence) integrate(constant float64) PolySequence {
	if len(ps) == 0 {
		// Shouldn't be relevant in practice, still treat it gracefully.
		return nil
	}

	result := make(PolySequence, len(ps))
	result[0] = ps[0].integrate(constant)

	for i := 1; i < len(ps); i++ {
		prev, err := result[i-1].Eval(result[i-1].XE)

		if err != nil {
			log.Printf("Bug: evaluate polynomial at xE must not fail: %v", err)
		}

		result[i] = ps[i].integrate(prev)
	}

	return result
}

// multiply multiplies all polynomial coefficients with factor in place. Factor must not be zero.
// Polynomial domains are unchanged.
func (ps PolySequence) multiply(factor float64) {
	for i := 0; i < len(ps); i++ {
		for j := 0; j < len(ps[i].Coeff); j++ {
			ps[i].Coeff[j] *= factor
		}
	}
}
