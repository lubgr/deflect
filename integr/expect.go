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
	Primary(d []deflect.NodalValue, t *testing.T)
	Reaction(r []deflect.NodalValue, t *testing.T)
	Failure(err error, t *testing.T)
	Interpolation(t *testing.T)
}

type noopExpectation struct{}

func (e *noopExpectation) Primary(d []deflect.NodalValue, t *testing.T)  {}
func (e *noopExpectation) Reaction(r []deflect.NodalValue, t *testing.T) {}
func (e *noopExpectation) Interpolation(t *testing.T)                    {}
func (e *noopExpectation) Failure(err error, t *testing.T)               {}

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

func (e *nodalExpectation) Reaction(r []deflect.NodalValue, t *testing.T) {
	t.Helper()
	e.nodal("reaction", r, e.reactions, e.tolReaction, t)
}

func (e *nodalExpectation) Primary(d []deflect.NodalValue, t *testing.T) {
	t.Helper()
	e.nodal("primary", d, e.primary, e.tolPrimary, t)
}

func (e *nodalExpectation) nodal(
	desc string,
	actual, expected []deflect.NodalValue,
	tol float64,
	t *testing.T,
) {
	t.Helper()

	for _, exp := range expected {
		idx := exp.Index
		pos := slices.IndexFunc(
			actual,
			func(result deflect.NodalValue) bool { return idx == result.Index },
		)

		if pos == -1 {
			t.Errorf("Nodal %v at %v/%v not found in %v results", desc, idx.NodalID, idx.Dof, len(actual))
			continue
		}

		approxEqual := scalar.EqualWithinAbsOrRel(actual[pos].Value, exp.Value, 1e-8, tol)

		if !approxEqual {
			t.Errorf(
				"Expected %v at %v/%v to be %v, got %v",
				desc,
				idx.NodalID,
				idx.Dof,
				exp.Value,
				actual[pos].Value,
			)
		}
	}
}

type nodalValues map[string]float64

type expectedInterpolationDescription struct {
	Element string
	Degree  uint
	// Limits the interpolation described to this range on the primary axis of the element. E.g., to
	// describe a constant solution in the first 1/3 of the elements lengths, Range is [0, 1/3].
	Range []float64
	// Specifies a sequence of x/value pairs, e.g. [[0, 5.0], [1.5 3]] would expected the
	// interpolation to yield 5 at x = 0 and 3 at x = 1.5.
	Values [][]float64
}

type expectedDescription struct {
	Tolerance     struct{ Primary, Reaction *float64 }
	Primary       map[string][]nodalValues
	Reaction      map[string][]nodalValues
	Interpolation expectedInterpolationDescription
	// A regular expression for the error description. If this field is not specified, success is
	// assumed and tested for.
	Failure *string
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

	nodal := newNodalExpectation(tmp.Expected.Tolerance.Primary, tmp.Expected.Tolerance.Reaction)
	var result []Expectation
	var errPrimary, errReact error

	if tmp.Expected.Failure == nil {
		result = append(result, &noFailureExpectation{})
	} else if pattern, err := regexp.Compile(*tmp.Expected.Failure); err != nil {
		return result, fmt.Errorf("couldn't compile '%v' to regexp", *tmp.Expected.Failure)
	} else {
		result = append(result, &regexFailureExpectation{pattern: pattern})
	}

	nodal.primary, errPrimary = translateToNodalExpectations(tmp.Expected.Primary)
	nodal.reactions, errReact = translateToNodalExpectations(tmp.Expected.Reaction)

	if err := errors.Join(errPrimary, errReact); err != nil {
		return nil, fmt.Errorf("failed to build nodal expectations: %w", err)
	} else if len(nodal.primary)+len(nodal.reactions) > 0 {
		result = append(result, &nodal)
	}

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
