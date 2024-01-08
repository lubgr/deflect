package xyz

import (
	"cmp"
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

// CompareIndices compares to Index instances in the spirit of [cmd.Compare].
func CompareIndices(i, j Index) int {
	if ids := cmp.Compare(i.NodalID, j.NodalID); ids != 0 {
		return ids
	}

	return cmp.Compare(string(i.Dof), string(j.Dof))
}
