package deflect

import (
	"maps"
	"testing"
)

func TestCompareIndices(t *testing.T) {
	cases := []struct {
		lhs      Index
		r        Index
		expected int
	}{
		{Index{"0", Ux}, Index{"5", Ux}, -1},
		{Index{"5", Ux}, Index{"5", Ux}, 0},
		{Index{"6", Ux}, Index{"5", Ux}, 1},
		{Index{"5", Ux}, Index{"5", Uz}, -1},
		{Index{"5", Phix}, Index{"5", Phiy}, 1},
		{Index{"abc", Phix}, Index{"xyz", Phiy}, -1},
	}

	for _, test := range cases {
		if actual := compareIndices(test.lhs, test.r); actual != test.expected {
			t.Errorf("Expected comparison between %v and %v to be %v, got %v",
				test.lhs, test.r, test.expected, actual)
		}
	}
}

func dofSet(data []Index) map[Index]struct{} {
	result := map[Index]struct{}{}

	for _, arg := range data {
		result[arg] = struct{}{}
	}

	return result
}

func TestIndexMapCreationWithoutOverlap(t *testing.T) {
	aux := Index{NodalID: "A", Dof: Ux}
	buz := Index{NodalID: "B", Dof: Uz}
	cphix := Index{NodalID: "C", Dof: Phix}
	cuz := Index{NodalID: "C", Dof: Uz}

	cases := [...]struct {
		free        []Index
		constrained []Index
		expected    map[Index]int
	}{
		{
			free:        []Index{},
			constrained: []Index{},
			expected:    map[Index]int{},
		},
		{
			free:        []Index{aux},
			constrained: []Index{},
			expected:    map[Index]int{aux: 0},
		},
		{
			free:        []Index{aux},
			constrained: []Index{aux},
			expected:    map[Index]int{aux: 0},
		},
		{
			free:        []Index{aux, buz, cphix, cuz},
			constrained: []Index{},
			expected:    map[Index]int{aux: 0, buz: 1, cuz: 2, cphix: 3},
		},
		{
			free:        []Index{},
			constrained: []Index{cuz, aux, cphix, cuz, buz},
			expected:    map[Index]int{aux: 0, buz: 1, cuz: 2, cphix: 3},
		},
		{
			free:        []Index{cuz, aux, cphix, cuz, buz},
			constrained: []Index{cphix, cuz},
			expected:    map[Index]int{cuz: 0, cphix: 1, aux: 2, buz: 3},
		},
		{
			free:        []Index{cuz, aux, cphix, cuz, buz},
			constrained: []Index{aux, cuz},
			expected:    map[Index]int{aux: 0, cuz: 1, buz: 2, cphix: 3},
		},
	}

	for _, test := range cases {
		free, constrained := dofSet(test.free), dofSet(test.constrained)
		actual := freeAndConstraintsToIndices(free, constrained)

		if !maps.Equal(actual, test.expected) {
			t.Errorf(
				"Given free/constrained indices %v/%v, expect mapping %v but got %v",
				test.free,
				test.constrained,
				test.expected,
				actual,
			)
		}
	}
}
