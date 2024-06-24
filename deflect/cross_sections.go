package deflect

import (
	"errors"
	"fmt"
)

type constantsCrossSection struct {
	area, iyy, izz, roll float64
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
	cs.izz = param["Izz"]

	for _, key := range []string{"A", "Iyy", "Izz"} {
		v, ok := param[key]
		if !ok {
			err = errors.Join(err, fmt.Errorf("cross section constant '%v' not found", key))
		} else if v <= 0 {
			err = errors.Join(err, fmt.Errorf("got non-positive cross section constant '%v' = %v", key, v))
		}
	}

	cs.roll = param["roll"] // Default zero if not given is desired

	return &cs, err
}

func (c *constantsCrossSection) Area() float64 {
	return c.area
}

func (c *constantsCrossSection) Iyy() float64 {
	return c.iyy
}

func (c *constantsCrossSection) Izz() float64 {
	return c.izz
}

func (c *constantsCrossSection) RollAngle() float64 {
	return c.roll
}

type rectangular struct {
	b, h, roll float64
}

// NewRectangularCrossSection instantiates a new rectangular cross section. Returns an error if the
// dimensions are not positive.
func NewRectangularCrossSection(b, h, roll float64) (CrossSection, error) {
	if b <= 0 || h <= 0 {
		return nil, fmt.Errorf("b/h must be positive, not %v/%v", b, h)
	} else if roll < 0 {
		return nil, fmt.Errorf("roll angle must not be negative, got %v", roll)
	}

	return &rectangular{b: b, h: h, roll: roll}, nil
}

func (r *rectangular) Area() float64 {
	return r.b * r.h
}

func (r *rectangular) Iyy() float64 {
	return r.b * r.h * r.h * r.h / 12.0
}

func (r *rectangular) Izz() float64 {
	return r.h * r.b * r.b * r.b / 12.0
}

func (r *rectangular) RollAngle() float64 {
	return r.roll
}
