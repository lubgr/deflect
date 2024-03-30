package deflect

import "testing"

func TestDofStringer(t *testing.T) {
	actual := Ux.String()
	if actual != "Ux" {
		t.Errorf("Expected Dof's String() to be 'Ux', got '%v'", actual)
	}
}

func TestDofSortOrder(t *testing.T) {
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
