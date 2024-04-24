package deflect

import (
	"slices"
	"strconv"
	"testing"
)

func TestTransformIntToString(t *testing.T) {
	cases := []struct {
		input  []int
		expect []string
	}{
		{input: []int{}, expect: []string{}},
		{input: []int{0}, expect: []string{"0"}},
		{input: []int{1, 2, 3}, expect: []string{"1", "2", "3"}},
	}

	for _, test := range cases {
		actual := transform(strconv.Itoa, test.input)

		if !slices.Equal(actual, test.expect) {
			t.Errorf("Failed int to string mapping, expected %v, got %v", test.expect, actual)
		}
	}
}
