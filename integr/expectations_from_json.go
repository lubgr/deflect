package integr

import (
	"encoding/json"
	"errors"
	"fmt"
	"regexp"
	"slices"

	"github.com/lubgr/deflect/deflect"
)

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

// ExpectationsFromJSON parses the given JSON data and constructs expectations that implement
// testing assertions from it.
func ExpectationsFromJSON(data []byte) ([]Expectation, error) {
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
