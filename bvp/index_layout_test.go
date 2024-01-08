package bvp

import (
	"maps"
	"testing"

	"github.com/lubgr/deflect/xyz"
)

func dofSet(data []xyz.Index) map[xyz.Index]struct{} {
	result := map[xyz.Index]struct{}{}

	for _, arg := range data {
		result[arg] = struct{}{}
	}

	return result
}

func TestIndexMapCreationWithoutOverlap(t *testing.T) {
	aux := xyz.Index{NodalID: "A", Dof: xyz.Ux}
	buz := xyz.Index{NodalID: "B", Dof: xyz.Uz}
	cphix := xyz.Index{NodalID: "C", Dof: xyz.Phix}
	cuz := xyz.Index{NodalID: "C", Dof: xyz.Uz}

	cases := [...]struct {
		free        []xyz.Index
		constrained []xyz.Index
		expected    map[xyz.Index]int
	}{
		{
			free:        []xyz.Index{},
			constrained: []xyz.Index{},
			expected:    map[xyz.Index]int{},
		},
		{
			free:        []xyz.Index{aux},
			constrained: []xyz.Index{},
			expected:    map[xyz.Index]int{aux: 0},
		},
		{
			free:        []xyz.Index{aux},
			constrained: []xyz.Index{aux},
			expected:    map[xyz.Index]int{aux: 0},
		},
		{
			free:        []xyz.Index{aux, buz, cphix, cuz},
			constrained: []xyz.Index{},
			expected:    map[xyz.Index]int{aux: 0, buz: 1, cuz: 2, cphix: 3},
		},
		{
			free:        []xyz.Index{},
			constrained: []xyz.Index{cuz, aux, cphix, cuz, buz},
			expected:    map[xyz.Index]int{aux: 0, buz: 1, cuz: 2, cphix: 3},
		},
		{
			free:        []xyz.Index{cuz, aux, cphix, cuz, buz},
			constrained: []xyz.Index{cphix, cuz},
			expected:    map[xyz.Index]int{cuz: 0, cphix: 1, aux: 2, buz: 3},
		},
		{
			free:        []xyz.Index{cuz, aux, cphix, cuz, buz},
			constrained: []xyz.Index{aux, cuz},
			expected:    map[xyz.Index]int{aux: 0, cuz: 1, buz: 2, cphix: 3},
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
