package xyz

import (
	"testing"
)

func TestStringer(t *testing.T) {
	actual := Ux.String()
	if actual != "u_x" {
		t.Errorf("Expected Dof's String() to be 'u_x', got '%v'", actual)
	}
}

func TestSortOrder(t *testing.T) {
	if !(Ux < Uz) {
		t.Errorf("Intuitive ordering Ux < Uz not satisfied")
	}
	if !(Uz < Phiy) {
		t.Errorf("Intuitive ordering Uz < Phiy not satisfied")
	}
	if !(Phiy < Phix) {
		t.Errorf("Intuitive ordering Phiy < Phix not satisfied")
	}
}

func TestCompareIndices(t *testing.T) {
	cases := []struct {
		lhs      Index
		rhs      Index
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
		if actual := CompareIndices(test.lhs, test.rhs); actual != test.expected {
			t.Errorf("Expected comparison between %v and %v to be %v, got %v",
				test.lhs, test.rhs, test.expected, actual)
		}
	}
}
