package deflect

// LinearElastic "implements" a linear-elastic material law by providing access to required material
// parameters. This does not encapsulate or abstract away the actual material law - an element
// formulation is expected to inline the material law when computing residual and tangent. This
// might change in the future when we add non-linear material laws.
type LinearElastic struct {
	// Young's modulus in N/m^2.
	YoungsModulus float64
	// Poisson's ratio, dimensionless. Only one, since we assume an isotropic material.
	PoissonsRatio float64
	// Density given in kg/m^3. We assume the Density is constant across every element.
	Density float64
}

func (el *LinearElastic) ShearModulus() float64 {
	return (el.YoungsModulus / (1.0 + el.PoissonsRatio)) / 2.0
}
