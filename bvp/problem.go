package bvp

import (
	"github.com/lubgr/deflect/elmt"
	"github.com/lubgr/deflect/xyz"
	"gonum.org/v1/gonum/mat"
)

// NodalValue associates a scalar solution value with an [xyz.Index].
type NodalValue struct {
	xyz.Index
	Value float64
}

// Problem stores data to solve a user specified boundary value problem.
type Problem struct {
	Nodes    []xyz.Node
	Elements []elmt.Element
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
		idx IndexLayout,
		strategy EquationSolver,
	) (primary, reactions []NodalValue, err error)
}
