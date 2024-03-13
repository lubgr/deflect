package deflect

import "fmt"

type neumannElementConcentrated struct {
	kind            Dof
	position, value float64
}

// NewElementConcentratedLoad instantiates an element Neumann boundary condition that applies a
// concentrated quantity at the given position, e.g. a force or torque.
func NewElementConcentratedLoad(kind Dof, pos, value float64) (NeumannElementBC, error) {
	if pos < 0 {
		return nil, fmt.Errorf("can't instantiate element force at negative x = %v", pos)
	}

	return &neumannElementConcentrated{kind: kind, position: pos, value: value}, nil
}

type neumannElementConstant struct {
	kind  Dof
	value float64
}

// NewElementConstantLoad instantiates an element Neumann boundary condition that applies a
// distributed, constant load.
func NewElementConstantLoad(kind Dof, value float64) NeumannElementBC {
	return &neumannElementConstant{kind: kind, value: value}
}

type neumannElementLinear struct {
	kind        Dof
	first, last float64
}

// NewElementLinearLoad instantiates an element Neumann boundary condition that applies a
// distributed, linear load.
func NewElementLinearLoad(kind Dof, first, last float64) NeumannElementBC {
	return &neumannElementLinear{kind: kind, first: first, last: last}
}
