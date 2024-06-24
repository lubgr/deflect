package deflect

import (
	"gonum.org/v1/gonum/mat"
)

// Dof describes a degree of freedom symbolically and with an intuitive ordering. The predefined
// constants shall be used.
//
//go:generate go run golang.org/x/tools/cmd/stringer -type=Dof
type Dof uint8

// Pre-defined degrees of freedom.
const (
	Ux Dof = iota
	Uz
	Uy
	Phiy
	Phiz
	Phix
)

// Fct denotes a function to be interpolated on the element level. These can be primary solution
// variables (nodal values) as well as stresses/internal forces.
//
//go:generate go run golang.org/x/tools/cmd/stringer -type=Fct -trimprefix Fct
type Fct uint8

// Pre-defined functions of interest that can be interpolated.
const (
	FctUnknown Fct = iota
	FctUx
	FctUz
	FctUy
	FctPhiy
	FctPhiz
	FctPhix
	FctNx
	FctVz
	FctVy
	FctMy
	FctMz
	FctMx
)

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
	Area() float64
	Iyy() float64
	Izz() float64
	RollAngle() float64
}

// NeumannElementBC is an opaque handle to be downcast by element implementations. It is always
// instantiated with a pointer.
type NeumannElementBC any

// Interpolation associates an element ID with a [Fct] descriptor for the quantity, and a variable
// number of piecewise polynomials. When there is only one polynomial, it spans the entire element
// domain (usually from zero to the length of the element). When there is more than one polynomial,
// there are sub-ranges/intervals. There will not be gaps in between these sub-ranges, and there
// will only be a single [PolyPiece] instance per sub-range.
type Interpolation struct {
	Element   string
	Quantity  Fct
	Piecewise PolySequence
}

// Element defines the common API of any finite element formulation implemented in this package.
type Element interface {
	// Assemble adds entries to the given tangent k and residual r, using the current primary nodal
	// values in d. Only non-linear elements need to read from primary. The matrices are global
	// entities, and the element uses indices to map symbolic indices to plain matrix indices, which
	// can be used to access the global matrices.
	Assemble(indices EqLayout, k *mat.SymDense, r, d *mat.VecDense)
	// AddLoad stores the given Neumann boundary condition to be used by the [Element.Assemble]
	// implementation. The same load can be added multiple times. The boolean return value indicates
	// if the load could be applied. If the inability to apply a BC is an error is up to the caller.
	AddLoad(bc NeumannElementBC) bool
	// RemoveLoad removes the given load. If there are multiple instances of the same load instance,
	// they are all removed. Removing a load that has previously not been added has no effect.
	RemoveLoad(bc NeumannElementBC)
	// Interpolate returns a sequence of piecewise polynomials to describe the primary solution or
	// stress/internal force in local element coordinates. Post-processing should be performed using
	// [PolySequence.TrimTrailingZeros] and [PolySequence.CompactIdentical] in that order, using a
	// context-dependent zero/equality tolerance (that's why it is not performed automatically). If
	// the element doesn't relate to the given dof, nil is returned.
	Interpolate(indices EqLayout, which Fct, d *mat.VecDense) PolySequence
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
	) (ProblemResult, error)
}

// ProblemResult is the API to retrieve BVP results as needed. An implementation can choose to
// evaluate individual results lazily or in batches, whatever makes most sense for the
// representation of the data.
type ProblemResult interface {
	Primary(i Index) (NodalValue, error)
	Reaction(i Index) (NodalValue, error)
	PrimaryAll() []NodalValue
	ReactionAll() []NodalValue
	Interpolate(elmtID string, quantity Fct, zeroTol float64) (Interpolation, error)
	InterpolateAll(zeroTol float64) []Interpolation
	Dimension() (total, net int)
}
