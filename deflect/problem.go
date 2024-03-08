package deflect

import (
	"gonum.org/v1/gonum/mat"
)

// Dof describes a degree of freedom symbolically and with a prefix to help maintain an intuitive
// ordering, e.g. horizontal displacement 'u_x' is listed before a rotation 'phi_y'. Dof instances
// should not be created by clients, use the predefined constants instead.
type Dof string

// Pre-defined degrees of freedom. The sort order prefixes should be treated as implementation
// details.
const (
	Ux   Dof = "00_u_x"
	Uz   Dof = "01_u_z"
	Uy   Dof = "02_u_y"
	Phiy Dof = "03_phi_y"
	Phiz Dof = "04_phi_z"
	Phix Dof = "05_phi_x"
)

// String returns the Dof description without the leading sort order prefix.
func (dof Dof) String() string {
	return string(dof)[3:]
}

// Node describes a mesh vertex by 3-dimensional coordinates and an identifier. Coordinates are
// always in meters.
type Node struct {
	ID string
	X  float64
	Y  float64
	Z  float64
}

// Index pairs a nodal id with a degree of freedom and is the complete symbolic representation of a
// quantity within an array of a boundary value problem. Dealing with raw integral matrix indices is
// error-prone, so Index is used as long as possible instead.
type Index struct {
	NodalID string
	Dof     Dof
}

// NodalValue associates a scalar solution value with an [Index].
type NodalValue struct {
	Index
	Value float64
}

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

// CrossSection provides access to all relevant geometric properties of a beam, truss, or frame
// cross section.
type CrossSection interface {
	// Area returns the cross section's area in m^2
	Area() float64
}

// Element defines the common API of any finite element formulation implemented in this package.
type Element interface {
	// Assemble adds entries to the given tangent k and residual r, using the current primary nodal
	// values in d. Only non-linear elements need to read from primary. The matrices are global
	// entities, and the element uses indices to map symbolic indices to plain matrix indices, which
	// can be used to access the global matrices.
	Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense)
	// Indices adds the indices this element uses (pairing of degree of freedom and nodal id) to set.
	Indices(set map[Index]struct{})
	// NumNodes returns the number of nodes this element is connected to.
	NumNodes() uint
	// ID returns the potentially user-facing id (hence it's a string).
	ID() string
}

// Transformer mutate the system of equation before and after it is solved. For example, an inclined
// support can be split into a Transformer and a Dirichlet BC, such that they act independently to
// link degrees of freedom through a trigonometric relation. Transformer instances don't prescribe
// values in r or d (this is done by NodalBC instances).
type Transformer interface {
	Pre(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense)
	Post(indices EqLayout, r, d *mat.VecDense)
}

// Problem stores data to solve a user specified boundary value problem.
type Problem struct {
	Nodes    []Node
	Elements []Element
	// Implementation note: nodal BCs are saved separately from basic.Node, since there is no coupling
	// between nodes and nodal boundary conditions (this is different from element loading, which is
	// tightly coupled).
	Neumann      []NodalValue // Only nodal Neumann BCs, no element loading
	Dirichlet    []NodalValue
	EqTransforms []Transformer
}

// EquationSolver implements an algorithm to solve a linear system of equations with a symmetric
// positive definite coefficient matrix. Examples: Cholesky or LU/QR decomposition, or an iterative
// method.
type EquationSolver interface {
	SolveLinearSystem(a mat.Symmetric, b, x *mat.VecDense) error
}

// ProblemSolver solves the given boundary value problem and returns the complete set of primary
// nodal values and reactions for nodes with Dirichlet BC.
type ProblemSolver interface {
	Solve(
		p *Problem,
		idx EqLayout,
		strategy EquationSolver,
	) (d, r []NodalValue, err error)
}
