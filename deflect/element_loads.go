package deflect

import "fmt"

type neumannConcentrated struct {
	kind            Dof
	position, value float64
}

// NewElementConcentratedLoad instantiates an element Neumann boundary condition that applies a
// concentrated quantity at the given position, e.g. a force or torque.
func NewElementConcentratedLoad(kind Dof, pos, value float64) (NeumannElementBC, error) {
	if pos < 0 {
		return nil, fmt.Errorf("can't instantiate element force at negative x = %v", pos)
	}

	return &neumannConcentrated{kind: kind, position: pos, value: value}, nil
}

type neumannConstant struct {
	kind  Dof
	value float64
}

// NewElementConstantLoad instantiates an element Neumann boundary condition that applies a
// distributed, constant load.
func NewElementConstantLoad(kind Dof, value float64) NeumannElementBC {
	return &neumannConstant{kind: kind, value: value}
}

type neumannLinear struct {
	kind        Dof
	first, last float64
}

// NewElementLinearLoad instantiates an element Neumann boundary condition that applies a
// distributed, linear load.
func NewElementLinearLoad(kind Dof, first, last float64) NeumannElementBC {
	return &neumannLinear{kind: kind, first: first, last: last}
}

// loadDispatch implements an exhaustive type switch over all NeumannElementBC types and calls the
// corresponding callback. This API shall be used instead of spreading identical type switches where
// needed. It gives us one single place to change when a new element load type is added/removed, and
// call sites can't be missed due to compiler errors. Therefore, we accept that this function
// signature is ugly and might even grow over time. The callbacks are expected to have side effects
// in their closure state.
//
// Implementation note: the following other options were considered and discarded.
// - New interface type with one method per element load type. This seems hideously verbose, in
// particular for scenarios with very little logic per load type.
// - Custom linter that inspects all type switches over NeumannElementBC. Seems more complicated
// than this approach.
// - Distinct struct with callback fields. This has no exhaustiveness guarantee at compile time,
// since it still allows default initialisation of callbacks, which are then unusable at runtime.
func loadDispatch(
	bc NeumannElementBC,
	concentrated func(*neumannConcentrated),
	constant func(*neumannConstant),
	linear func(*neumannLinear),
) {
	switch load := bc.(type) {
	case *neumannConcentrated:
		concentrated(load)
	case *neumannConstant:
		constant(load)
	case *neumannLinear:
		linear(load)
	}
}
