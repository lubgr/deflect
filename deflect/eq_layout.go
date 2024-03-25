package deflect

import (
	"cmp"
	"fmt"
	"slices"
	"strings"
)

// EqLayout manages the mapping between symbolic and integral matrix indices. This is used for
// elements and boundary conditions to know where in the global matrices/vectors they should add
// their coefficients during assembly. It collects lookup errors behind the scene. These errors must
// be manually checked at a reasonable point in time - the individual lookups will _not_ fail, but
// return zero indices instead. This means an out-of-bounds panic is avoided, but any matrices
// populated with the erroneous index mapping can't be used. In the context of tangent/residual
// assembly, the solver process should probably be aborted.
type EqLayout struct {
	indices  map[Index]int
	inverse  []Index
	failures int
	failed   []string
}

// eqSize returns the total size of the system of equations, including Dirichlet-constrained nodes.
func (il *EqLayout) eqSize() int {
	return len(il.indices)
}

// mapOne looks up and returns the mapped-to plain matrix index. Returns zero if isym is not mapped
// and internally accumulates an error, see [EqLayout.failure] and [EqLayout].
func (il *EqLayout) mapOne(isym Index) (i int) {
	var ok bool
	i, ok = il.indices[isym]
	il.saveFailure(ok, isym)
	return i
}

// mapTwo is identical to [mapOne], but with two lookups.
func (il *EqLayout) mapTwo(isym, jsym Index) (i, j int) {
	return il.mapOne(isym), il.mapOne(jsym)
}

// mapThree is identical to [mapOne], but with three lookups.
func (il *EqLayout) mapThree(isym, jsym, ksym Index) (i, j, k int) {
	return il.mapOne(isym), il.mapOne(jsym), il.mapOne(ksym)
}

// mapFour is identical to [mapOne], but with four lookups.
func (il *EqLayout) mapFour(isym, jsym, ksym, lsym Index) (i, j, k, l int) {
	return il.mapOne(isym), il.mapOne(jsym), il.mapOne(ksym), il.mapOne(lsym)
}

// unmap looks up the reverse mapping, but does not implement any error handling - i must have been
// retrieved by lookup functions like [mapOne] beforehand.
func (il *EqLayout) unmap(i int) (isym Index) {
	return il.inverse[i]
}

func (il *EqLayout) saveFailure(ok bool, symbolic Index) {
	if !ok {
		il.failures++

		if len(il.failed) < 5 {
			il.failed = append(il.failed, fmt.Sprintf("%v/%v", symbolic.NodalID, symbolic.Dof))
		}
	}
}

// failure returns a non-nil error if any of the preceding lookup operations failed. The internal
// error state is _not_ flushed. In order to clear the error, a new EqLayout instance must be
// instantiated.
func (il *EqLayout) failure() error {
	if il.failures == 0 {
		return nil
	}

	list := strings.Join(il.failed, ", ")

	if il.failures > len(il.failed) {
		list += ", ..."
	}

	return fmt.Errorf("%v index lookup failure(s) (%v)", il.failures, list)

}

// NewEqLayout creates a new index layout for the given boundary value problem.
func NewEqLayout(p *Problem) (EqLayout, error) {
	indices, err := createIndexMap(p.Elements, p.Dirichlet)

	return newEqLayoutDirect(indices), err
}

func newEqLayoutDirect(indices map[Index]int) EqLayout {
	return EqLayout{indices: indices, inverse: invertIndexMap(indices)}
}

// IndexMap constructs the mapping from symbolic to plain matrix indices. It returns an error when
// there are duplicate Dirichlet indices.
func createIndexMap(elements []Element, dirichlet []NodalValue) (map[Index]int, error) {
	indices := map[Index]struct{}{}
	constrained := map[Index]struct{}{}
	var duplicates []string

	for _, d := range dirichlet {
		if _, ok := constrained[d.Index]; ok {
			duplicates = append(duplicates, fmt.Sprintf("%v/%v", d.NodalID, d.Dof))
		}
		constrained[d.Index] = struct{}{}
	}

	for _, element := range elements {
		element.Indices(indices)
	}

	var err error

	if len(duplicates) > 0 {
		err = fmt.Errorf("duplicate Dirichlet BCs %v", strings.Join(duplicates, ", "))
	}

	return freeAndConstraintsToIndices(indices, constrained), err
}

// freeAndConstraintsToIndices constructs a lookup table to map symbolic indices to plain array
// indices. Constrained indices override free indices. The mapping orders indices such that matrices
// and vectors that are populated with this mapping are partitioned according to the given
// constraints. Constrained symbolic indices have plain top/left matrix indices, unconstrained ones
// right/bottom.
func freeAndConstraintsToIndices(
	free map[Index]struct{},
	constrained map[Index]struct{},
) map[Index]int {
	indices := make([]Index, 0, len(free))

	for index := range constrained {
		indices = append(indices, index)
	}

	for index := range free {
		if _, ok := constrained[index]; !ok {
			indices = append(indices, index)
		}
	}

	part := len(constrained)
	slices.SortFunc(indices[0:part], compareIndices)
	slices.SortFunc(indices[part:], compareIndices)

	result := make(map[Index]int)

	for i, index := range indices {
		result[index] = i
	}

	return result
}

func invertIndexMap(from map[Index]int) []Index {
	to := make([]Index, len(from))

	for symbolic, plain := range from {
		to[plain] = symbolic
	}

	return to
}

// compareIndices compares to Index instances in the spirit of [cmd.Compare].
func compareIndices(i, j Index) int {
	if ids := cmp.Compare(i.NodalID, j.NodalID); ids != 0 {
		return ids
	}

	return cmp.Compare(string(i.Dof), string(j.Dof))
}
