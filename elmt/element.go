package elmt

import (
	"github.com/lubgr/deflect/xyz"
	"gonum.org/v1/gonum/mat"
)

// Element defines the common API of any finite element formulation implemented in this package.
type Element interface {
	// Assemble adds entries to the given tangent k and residual rhs, using the current primary nodal
	// values in primary. Only non-linear elements need to read from primary. The matrices are global
	// entities, and the element uses indices to map symbolic indices to plain matrix indices, which
	// can be used to access the global matrices.
	Assemble(indices map[xyz.Index]int, k *mat.SymDense, rhs, primary *mat.VecDense) error
	// Indices adds the indices this element uses (pairing of degree of freedom and nodal id) to set.
	Indices(set map[xyz.Index]struct{})
	// NumNodes returns the number of nodes this element is connected to.
	NumNodes() uint
	// ID returns the potentially user-facing id (hence it's a string).
	ID() string
}
