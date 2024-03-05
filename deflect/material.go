package deflect

// Material defines the API to retrieve material parameters to be used in element formulations. We
// currently restrain ourselves to linear-elastic materials, and the constant material parameters
// are directly accessed by each element to compute its residual and tangent. The Material
// abstraction is hence intentionally thin for now, and it is fine that element formulations and
// material law are tightly coupled.
//
// In the future, non-linear material laws can be an important feature. Then, the Material API must
// be reworked (e.g. return a function to compute residual/tangent given current strain), and we
// would probably try to decouple the element's interpolation functionality from the material law.
// Then, an element formulation can work with different injected material instances.
type Material struct {
	CrossSection
	LinearElastic
}
