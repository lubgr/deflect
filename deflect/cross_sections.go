package deflect

import (
	"errors"
	"fmt"
	"math"
)

type constantsCrossSection struct {
	area, iyy, izz, ixx, roll float64
}

// NewConstantsCrossSections instantiates a cross section with all parameters specified as
// constants. The following keys are expected:
// - A
// - Iyy
// - Izz
// Optional parameters are
// - Ixx
// - roll (angle)
// Returns an error if any of the mandatory parameters is not positive or can't be found in param.
// An error is also returned if an optional parameter is negative. No error is returned if an
// optional parameter is zero.
func NewConstantsCrossSections(param map[string]float64) (CrossSection, error) {
	var cs constantsCrossSection
	var err error

	cs.area = param["A"]
	cs.iyy = param["Iyy"]
	cs.izz = param["Izz"]

	for _, key := range []string{"A", "Iyy", "Izz"} {
		v, ok := param[key]
		if !ok {
			err = errors.Join(err, fmt.Errorf("cross section constant %v not found", key))
		} else if v <= 0 {
			err = errors.Join(err, fmt.Errorf("non-positive cross section constant %v = %v", key, v))
		}
	}

	// Default zeros if not given are desired
	cs.ixx = param["Ixx"]
	cs.roll = param["roll"]

	if cs.ixx < 0 {
		err = errors.Join(err, fmt.Errorf("negative cross section constant Ixx = %v", cs.ixx))
	}

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

func (c *constantsCrossSection) Ixx() float64 {
	return c.ixx
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

func (r *rectangular) Ixx() float64 {
	a, b := r.b, r.h

	// Make sure a denotes the longer side of the rectangle, while b is the shorter side
	if a < b {
		a, b = b, a
	}

	// Standard formula for rectangles, e.g. "Roark's Formulas for Stress and Strain" (W.C. Young,
	// R.G. Budynas), 7th ed., Table 10.1, p. 401: a*b³/16*(16/3 - 3.36*b/a*(1 - b⁴/(12*a⁴)))
	return a * math.Pow(b, 3) * (1.0/3.0 - 3.36/16.0*b/a*(1.0-math.Pow(b, 4)/(12.0*math.Pow(a, 4))))
}

func (r *rectangular) RollAngle() float64 {
	return r.roll
}
