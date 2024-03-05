package deflect

import (
	"cmp"
	"fmt"
	"slices"
	"strings"
)

// IndexLayout stores the mapping between symbolic and integral matrix indices. This is used for
// elements and boundary conditions to know where in the global matrices/vectors they should add
// their coefficients during assembly.
type IndexLayout struct {
	indices map[Index]int
	inverse []Index
}

// NewIndexLayout creates a new index layout for the given boundary value problem.
func NewIndexLayout(p *Problem) (IndexLayout, error) {
	var result IndexLayout
	var err error

	result.indices, err = createIndexMap(p.Elements, p.Dirichlet)
	result.inverse = invertIndexMap(result.indices)

	return result, err
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