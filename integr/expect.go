package integr

import (
	"encoding/json"
	"errors"
	"fmt"
	"regexp"
	"slices"
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

type expectedPolynomial struct {
	quantity deflect.Fct
	degree   int
	fctRange []float64
	boundary []float64
	eval     [][]float64
}

type interpolationExpectation struct {
	noopExpectation
	tolerance float64
	functions map[string][]expectedPolynomial
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

type expectedDescription struct {
	Tolerance     struct{ Primary, Reaction, Polynomial *float64 }
	Primary       map[string][]nodalValues
	Reaction      map[string][]nodalValues
	Interpolation map[string][]expectedInterpolationDescription
	// A regular expression for the error description. If this field is not specified, success is
	// assumed and tested for.
	Failure *string
}

type nodalValues map[string]float64

type expectedInterpolationDescription struct {
	Kind   string
	Degree int
	// If not nil or empty, limits the interpolation described to this range on the
	// primary axis of the element. E.g., to describe a constant solution in the first 1/3
	// of the elements lengths, Range is [0, 1/3]. If nil/empty, the range spans the
	// entire element length.
	Range []float64
	// When empty, nothing is tested. Length 1 makes sense for constant polynomials. Length 2 is for
	// all other polynomials of higher degree (linear onwards), the first value is compared against
	// the evaluation of the smallest x of the left-most piecewise polynomial, the second value
	// against the largest x of the right-most polynomial.
	Boundary []float64
	// Specifies a sequence of x/value pairs, e.g. [[0, 5.0], [1.5 3]] would expected the
	// interpolation to yield 5 at x = 0 and 3 at x = 1.5.
	Eval [][]float64
}

// ExpectFromJSON parses the given JSON data and constructs expectations that implement testing
// assertions from it.
func ExpectFromJSON(data []byte) ([]Expectation, error) {
	var tmp struct {
		Expected expectedDescription
	}

	if err := json.Unmarshal(data, &tmp); err != nil {
		return nil, fmt.Errorf("couldn't build expectations from JSON: %w", err)
	}

	var result []Expectation

	if tmp.Expected.Failure == nil {
		result = append(result, &noFailureExpectation{})
	} else if pattern, err := regexp.Compile(*tmp.Expected.Failure); err != nil {
		return result, fmt.Errorf("couldn't compile '%v' to regexp", *tmp.Expected.Failure)
	} else {
		result = append(result, &regexFailureExpectation{pattern: pattern})
	}

	nodal := newNodalExpectation(tmp.Expected.Tolerance.Primary, tmp.Expected.Tolerance.Reaction)
	var errPrimary, errReact error
	nodal.primary, errPrimary = translateToNodalExpectations(tmp.Expected.Primary)
	nodal.reactions, errReact = translateToNodalExpectations(tmp.Expected.Reaction)

	if err := errors.Join(errPrimary, errReact); err != nil {
		return nil, fmt.Errorf("failed to build nodal expectations: %w", err)
	} else if len(nodal.primary)+len(nodal.reactions) > 0 {
		result = append(result, &nodal)
	}

	interpolation, err := translateInterpolations(
		tmp.Expected.Interpolation,
		tmp.Expected.Tolerance.Polynomial,
	)

	if err != nil {
		return nil, fmt.Errorf("failed to build test expectations for interpolations: %w", err)
	}

	result = append(result, interpolation)

	return result, nil
}

func translateToNodalExpectations(from map[string][]nodalValues) ([]deflect.NodalValue, error) {
	dofs := map[string]deflect.Dof{
		"Fx":   deflect.Ux,
		"Fz":   deflect.Uz,
		"Fy":   deflect.Uy,
		"My":   deflect.Phiy,
		"Mz":   deflect.Phiz,
		"Mx":   deflect.Phix,
		"Ux":   deflect.Ux,
		"Uz":   deflect.Uz,
		"Uy":   deflect.Uy,
		"Phiy": deflect.Phiy,
		"Phiz": deflect.Phiz,
		"Phix": deflect.Phix,
	}
	expect := make([]deflect.NodalValue, 0, 2*len(from)) // Capacity is just a guess
	var err error

	for nodalID, desc := range from {
		for _, reaction := range desc {
			for name, value := range reaction {
				dof, ok := dofs[name]

				if !ok {
					err = errors.Join(err, fmt.Errorf("unknown reaction/degree of freedom '%v'", name))
					continue
				}

				idx := deflect.Index{NodalID: nodalID, Dof: dof}
				expect = append(expect, deflect.NodalValue{Index: idx, Value: value})
			}
		}
	}

	return expect, err
}

func translateInterpolations(
	from map[string][]expectedInterpolationDescription,
	tol *float64,
) (Expectation, error) {
	expect := interpolationExpectation{tolerance: 1e-8, functions: map[string][]expectedPolynomial{}}

	if tol != nil {
		expect.tolerance = *tol
	}

	for elmtID, descriptions := range from {
		perElement := expect.functions[elmtID]

		for _, desc := range descriptions {
			p, err := translateInterpolation(&desc)

			if err != nil {
				return nil, fmt.Errorf("construct interpolation for element %v: %w", elmtID, err)
			}

			perElement = append(perElement, p)
		}

		expect.functions[elmtID] = perElement
	}

	return &expect, nil
}

func translateInterpolation(from *expectedInterpolationDescription) (expectedPolynomial, error) {
	quantities := map[string]deflect.Fct{
		"Nx":   deflect.FctNx,
		"Vz":   deflect.FctVz,
		"Vy":   deflect.FctVy,
		"My":   deflect.FctMy,
		"Mz":   deflect.FctMz,
		"Mx":   deflect.FctMx,
		"Ux":   deflect.FctUx,
		"Uz":   deflect.FctUz,
		"Uy":   deflect.FctUy,
		"Phiy": deflect.FctPhiy,
		"Phiz": deflect.FctPhiz,
		"Phix": deflect.FctPhix,
	}

	var p expectedPolynomial

	q, ok := quantities[from.Kind]

	if !ok {
		return p, fmt.Errorf("unknown interpolation '%v'", from.Kind)
	}

	p.quantity = q
	p.degree = from.Degree

	if r := len(from.Range); r != 0 && r != 2 {
		return p, fmt.Errorf("interpolation range with %v values, must be 0 or 2", r)
	}

	p.fctRange = from.Range

	if slices.ContainsFunc(from.Eval, func(pair []float64) bool { return len(pair) != 2 }) {
		return p, errors.New("interpolation sample values must be pairs")
	}

	p.eval = from.Eval

	if len(from.Boundary) > 2 {
		return p, errors.New("at most 2 boundary values (left/right) can be specified")
	}

	p.boundary = from.Boundary

	return p, nil
}
