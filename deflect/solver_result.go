package deflect

import (
	"errors"
	"fmt"
	"log"
	"slices"

	"gonum.org/v1/gonum/mat"
)

type solverResult struct {
	total, net int
	d, r       *mat.VecDense
	indices    EqLayout
	elements   []Element
	// Group solution/reaction with symbolic index:
	dIndexed, rIndexed []NodalValue
}

func (sr *solverResult) Primary(i Index) (NodalValue, error) {
	plain := sr.indices.mapOne(i)
	return NodalValue{Index: i, Value: sr.d.AtVec(plain)}, sr.indices.flushFailure()
}

func (sr *solverResult) Reaction(i Index) (NodalValue, error) {
	plain := sr.indices.mapOne(i)
	return NodalValue{Index: i, Value: sr.r.AtVec(plain)}, sr.indices.flushFailure()
}

func (sr *solverResult) PrimaryAll() []NodalValue {
	sr.dIndexed = sr.indexPaired(sr.dIndexed, sr.d)
	return sr.dIndexed
}

func (sr *solverResult) ReactionAll() []NodalValue {
	sr.rIndexed = sr.indexPaired(sr.rIndexed, sr.r)
	return sr.rIndexed
}

func (sr *solverResult) indexPaired(paired []NodalValue, raw *mat.VecDense) []NodalValue {
	if paired != nil {
		return paired
	}

	paired = make([]NodalValue, sr.total)

	for i := 0; i < sr.total; i++ {
		symbolic := sr.indices.unmap(i)
		paired[i] = NodalValue{Index: symbolic, Value: raw.AtVec(i)}
	}

	if err := sr.indices.flushFailure(); err != nil {
		log.Panicf("This is a bug, unmapping 0...%v to symbolic indices failed: %v", sr.total, err)
	}

	return paired
}

func (sr *solverResult) Interpolate(
	elmtID string,
	quantity Fct,
	zeroTol float64,
) (Interpolation, error) {
	idx := slices.IndexFunc(sr.elements, func(e Element) bool { return e.ID() == elmtID })

	if idx == -1 {
		return Interpolation{}, fmt.Errorf("no element with ID '%v' found", elmtID)
	}

	interpolation := sr.elements[idx].Interpolate(sr.indices, quantity, sr.d)
	interpolation.TrimTrailingZeros(zeroTol)
	interpolation = interpolation.CompactIdentical(zeroTol)

	result := Interpolation{
		Element:   sr.elements[idx].ID(),
		Quantity:  quantity,
		Piecewise: interpolation,
	}

	return result, sr.indices.flushFailure()
}

func (sr *solverResult) InterpolateAll(zeroTol float64) []Interpolation {
	quantities := [...]Fct{
		FctUx,
		FctUz,
		FctUy,
		FctPhiy,
		FctPhiz,
		FctPhix,
		FctNx,
		FctVz,
		FctVy,
		FctMy,
		FctMz,
		FctMx,
	}
	result := make([]Interpolation, 0, len(quantities)*len(sr.elements))
	var err error

	for _, q := range quantities {
		for _, e := range sr.elements {
			poly, errSingle := sr.Interpolate(e.ID(), q, zeroTol)
			err = errors.Join(err, errSingle)

			if errSingle == nil && poly.Piecewise != nil {
				result = append(result, poly)
			}
		}
	}

	if err != nil {
		log.Printf("Bug: interpolating all polynomials must not fail: %v", err)
	}

	return result
}

func (sr *solverResult) Dimension() (total, net int) {
	return sr.total, sr.net
}
