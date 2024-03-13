package deflect

import "testing"

func TestConcentratedLoadConstruction(t *testing.T) {
	_, err := NewElementConcentratedLoad(Ux, -0.5, 123.0)

	if err == nil {
		t.Errorf("BC construction succeeded despite given a negative position")
	}

	_, err = NewElementConcentratedLoad(Ux, 0.5, -123.0)

	if err != nil {
		t.Errorf("BC construction failed despite given correct parameters")
	}
}
