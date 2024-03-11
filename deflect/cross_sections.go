package deflect

import (
	"errors"
	"fmt"
)

type constantsCrossSection struct {
	area, iyy float64
}

// NewConstantsCrossSections instantiates a cross section with all parameters specified as
// constants. The following keys are expected:
// - A
// - Iyy
// Returns an error if anything is not positive, or if the required parameters can't be found in
// param.
func NewConstantsCrossSections(param map[string]float64) (CrossSection, error) {
	var cs constantsCrossSection
	var err error

	cs.area = param["A"]
	cs.iyy = param["Iyy"]

	for _, key := range []string{"A", "Iyy"} {
		v, ok := param[key]
		if !ok {
			err = errors.Join(err, fmt.Errorf("cross section constant '%v' not found", key))
		} else if v <= 0 {
			err = errors.Join(err, fmt.Errorf("got non-positive cross section constant '%v' = %v", key, v))
		}
	}

	return &cs, err
}

func (c *constantsCrossSection) Area() float64 {
	return c.area
}

func (c *constantsCrossSection) Iyy() float64 {
	return c.iyy
}

type rectangular struct {
	b, h float64
}

// NewRectangularCrossSection instantiates a new rectangular cross section. Returns an error if the
// dimensions are not positive.
func NewRectangularCrossSection(b, h float64) (CrossSection, error) {
	if b <= 0.0 || h <= 0.0 {
		return nil, fmt.Errorf("b/h must be positive, not %v/%v", b, h)
	}

	return &rectangular{b: b, h: h}, nil
}

func (r *rectangular) Area() float64 {
	return r.b * r.h
}

func (r *rectangular) Iyy() float64 {
	return r.b * r.h * r.h * r.h / 12.0
}
