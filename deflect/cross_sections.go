package deflect

import "fmt"

type rectangular struct {
	b, h float64
}

func (r *rectangular) Area() float64 {
	return r.b * r.h
}

// NewRectangularCrossSection instantiates a new rectangular cross section. Returns an error if the
// dimensions are not positive.
func NewRectangularCrossSection(b, h float64) (CrossSection, error) {
	if b <= 0.0 || h <= 0.0 {
		return nil, fmt.Errorf("b/h must be positive, not %v/%v", b, h)
	}

	return &rectangular{b: b, h: h}, nil
}
