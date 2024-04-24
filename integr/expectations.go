package integr

import (
	"fmt"
	"regexp"
	"strings"
	"testing"

	"github.com/lubgr/deflect/deflect"
	"gonum.org/v1/gonum/floats/scalar"
)

// Expectation offers an interface for integration tests to run assertions after solving a boundary
// value problem.
type Expectation interface {
	Primary(r deflect.ProblemResult, t *testing.T)
	Reaction(r deflect.ProblemResult, t *testing.T)
	Interpolated(r deflect.ProblemResult, t *testing.T)
	Failure(err error, t *testing.T)
}

type noopExpectation struct{}

func (e *noopExpectation) Primary(r deflect.ProblemResult, t *testing.T)      {}
func (e *noopExpectation) Reaction(r deflect.ProblemResult, t *testing.T)     {}
func (e *noopExpectation) Interpolated(r deflect.ProblemResult, t *testing.T) {}
func (e *noopExpectation) Failure(err error, t *testing.T)                    {}

type noFailureExpectation struct {
	noopExpectation
}

func (e *noFailureExpectation) Failure(err error, t *testing.T) {
	t.Helper()
	if err != nil {
		t.Fatalf("Unexpected BVP failure: %v", err)
	}
}

type regexFailureExpectation struct {
	noopExpectation
	pattern *regexp.Regexp
}

func (e *regexFailureExpectation) Failure(err error, t *testing.T) {
	t.Helper()
	if err == nil {
		t.Errorf("Expected failure with message matching '%v', but got no error", e.pattern.String())
		return
	}

	msg := err.Error()
	lower := strings.ToLower(msg)

	if !(e.pattern.MatchString(msg) || e.pattern.MatchString(lower)) {
		t.Errorf("Error message '%v' didn't match expected '%v'", msg, e.pattern.String())
	}
}

type nodalExpectation struct {
	noopExpectation
	primary, reactions      []deflect.NodalValue
	tolPrimary, tolReaction float64
}

func newNodalExpectation(primary, reaction *float64) nodalExpectation {
	result := nodalExpectation{tolPrimary: 1e-8, tolReaction: 1e-8}

	if primary != nil {
		result.tolPrimary = *primary
	}
	if reaction != nil {
		result.tolReaction = *reaction
	}

	return result
}

func (e *nodalExpectation) Reaction(r deflect.ProblemResult, t *testing.T) {
	t.Helper()
	e.nodal("reaction", r.Reaction, e.reactions, e.tolReaction, t)
}

func (e *nodalExpectation) Primary(r deflect.ProblemResult, t *testing.T) {
	t.Helper()
	e.nodal("primary", r.Primary, e.primary, e.tolPrimary, t)
}

func (e *nodalExpectation) nodal(
	desc string,
	query func(deflect.Index) (deflect.NodalValue, error),
	expected []deflect.NodalValue,
	tol float64,
	t *testing.T,
) {
	t.Helper()

	for _, exp := range expected {
		idx := exp.Index
		actual, err := query(idx)

		if err != nil {
			t.Errorf("Couldn't lookup nodal %v at %v/%v: %v", desc, idx.NodalID, idx.Dof, err)
			continue
		}

		approxEqual := scalar.EqualWithinAbsOrRel(actual.Value, exp.Value, 1e-8, tol)

		if !approxEqual {
			t.Errorf(
				"Expected %v at %v/%v to be %v, got %v",
				desc,
				idx.NodalID,
				idx.Dof,
				exp.Value,
				actual.Value,
			)
		}
	}
}

type interpolationExpectation struct {
	noopExpectation
	tolerance float64
	functions map[string][]expectedPolynomial
}

type expectedPolynomial struct {
	quantity deflect.Fct
	degree   int
	fctRange []float64
	boundary []float64
	eval     [][]float64
}

func (e *interpolationExpectation) Interpolated(r deflect.ProblemResult, t *testing.T) {
	for elmtID, interpolations := range e.functions {
		for _, expect := range interpolations {
			p, err := e.findPolySegment(r, elmtID, &expect)

			if err != nil {
				t.Errorf("Couldn't find polynomial segment for %v/%v: %v", elmtID, expect.quantity, err)
				continue
			}

			e.testInterpolation(elmtID, &expect, p, t)
		}
	}
}

func (e *interpolationExpectation) findPolySegment(
	r deflect.ProblemResult,
	elmtID string,
	expect *expectedPolynomial,
) (deflect.PolyPiece, error) {
	poly, err := r.Interpolate(elmtID, expect.quantity, e.tolerance)

	if err != nil {
		return deflect.PolyPiece{}, fmt.Errorf("query interpolation: %w", err)
	}

	ps := poly.Piecewise

	// If no range has been given, the polynomial characteristics are expected to hold for the
	// entire element. There must not be any sub-domains, and hence only a single coefficient set.
	if expect.fctRange == nil && len(ps) == 1 {
		return ps[0], nil
	} else if expect.fctRange == nil && len(ps) != 1 {
		return deflect.PolyPiece{}, fmt.Errorf("%v polynomials on element %v", len(ps), elmtID)
	}

	// An explicit range has been specified, so we search for the corresponding interval of the
	// polynomial. The length of the range slice must be 2, this has been checked earlier.
	x0, xE := expect.fctRange[0], expect.fctRange[1]
	idx := ps.IntervalIndex(x0, xE)

	if idx == -1 {
		err := fmt.Errorf("couldn't find [%v, %v] in %v polynomial intervals", x0, xE, len(ps))
		return deflect.PolyPiece{}, err
	}

	return ps[idx], nil
}

func (e *interpolationExpectation) testInterpolation(
	elmtID string,
	expect *expectedPolynomial,
	p deflect.PolyPiece,
	t *testing.T,
) {
	what := fmt.Sprintf("%v/%v", elmtID, expect.quantity)
	deg := p.Degree()

	if deg != expect.degree {
		t.Errorf(
			"Expected %v polynomial degree to be %v, got %v: %v",
			what,
			expect.degree,
			deg,
			p.Coeff,
		)
	}

	if len(p.Coeff) == 0 {
		t.Errorf("Invalid polynomial with empty coefficients: %v", what)
		return
	}

	// The easiest way to test boundary values is to append them to other (probably often empty)
	// values to evaluate.
	if len(expect.boundary) == 1 {
		expect.eval = append(expect.eval, []float64{p.X0, expect.boundary[0]})
	} else if len(expect.boundary) == 2 {
		expect.eval = append(expect.eval,
			[]float64{p.X0, expect.boundary[0]}, []float64{p.XE, expect.boundary[1]})
	}

	for _, sample := range expect.eval {
		x, value := sample[0], sample[1]
		actual, err := p.Eval(x)

		if err != nil {
			t.Errorf("Couldn't evaluate polynomial: %v", err)
			continue
		}

		if !scalar.EqualWithinAbsOrRel(value, actual, e.tolerance, e.tolerance) {
			t.Errorf("Expected %v [%v, %v] at x=%v to be %v, got %v", what, p.X0, p.XE, x, value, actual)
		}
	}
}
