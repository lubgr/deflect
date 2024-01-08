package elmt

import (
	"errors"
	"fmt"

	"github.com/lubgr/deflect/xyz"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/spatial/r3"
)

// Truss2dDisjoint specifies how a 2d truss is connected to its first and second node.
type Truss2dDisjoint uint8

// Enumerators for 2d truss joints.
const (
	Truss2dDisjointNone Truss2dDisjoint = iota
	Truss2dDisjointFirstNode
	Truss2dDisjointSecondNode
)

type truss2d struct {
	id       string
	n0, n1   *xyz.Node
	material *Material
	disjoint Truss2dDisjoint
}

func (t *truss2d) Assemble(indices map[xyz.Index]int, k *mat.SymDense,
	rhs, primary *mat.VecDense) error {
	local := t.localTangent()
	rot := t.rotation()

	// Computes the global tangent k = rot^T * (rot^T * local)^T
	global := mat.NewDense(4, 4, nil)
	global.Mul(rot.T(), local)
	global.Mul(rot.T(), global.T())

	ux0, uz0, ux1, uz1, err := t.mappedIndices(indices)

	if err != nil {
		return fmt.Errorf("element tangent assembly: %w", err)
	}

	add := func(i, j int, value float64) {
		k.SetSym(i, j, k.At(i, j)+value)
	}

	add(ux0, ux0, global.At(0, 0))
	add(ux0, uz0, global.At(0, 1))
	add(ux0, ux1, global.At(0, 2))
	add(ux0, uz1, global.At(0, 3))

	add(uz0, uz0, global.At(1, 1))
	add(uz0, ux1, global.At(1, 2))
	add(uz0, uz1, global.At(1, 3))

	add(ux1, ux1, global.At(2, 2))
	add(ux1, uz1, global.At(2, 3))

	add(uz1, uz1, global.At(3, 3))

	return nil
}

func (t *truss2d) mappedIndices(lookup map[xyz.Index]int) (int, int, int, int, error) {
	mapWithErr := func(symbolic xyz.Index) (int, error) {
		plain, ok := lookup[symbolic]

		if !ok {
			return -1, fmt.Errorf("%v not found", symbolic)
		}

		return plain, nil
	}

	indices := t.indicesAsArray()
	ux0, err0 := mapWithErr(indices[0])
	uz0, err1 := mapWithErr(indices[1])
	ux1, err2 := mapWithErr(indices[2])
	uz1, err3 := mapWithErr(indices[3])

	return ux0, uz0, ux1, uz1, errors.Join(err0, err1, err2, err3)
}

func (t *truss2d) indicesAsArray() *[4]xyz.Index {
	indices := [...]xyz.Index{
		{NodalID: t.n0.ID, Dof: xyz.Ux},
		{NodalID: t.n0.ID, Dof: xyz.Uz},
		{NodalID: t.n1.ID, Dof: xyz.Ux},
		{NodalID: t.n1.ID, Dof: xyz.Uz},
	}

	return &indices
}

func (t *truss2d) localTangent() *mat.SymDense {
	k := mat.NewSymDense(4, nil)

	if t.disjoint != Truss2dDisjointNone {
		return k
	}

	l := length(t.n0, t.n1)
	EA := t.material.YoungsModulus * t.material.Area()

	k.SetSym(0, 0, EA/l)
	k.SetSym(0, 2, -EA/l)
	k.SetSym(2, 2, EA/l)

	return k
}

func (t *truss2d) rotation() *mat.Dense {
	// Gonum's spatial/r2 package doesn't seem to have a way to create a matrix representation of a 2d
	// rotation (only for 3d), so we do it manually.
	line := r3.Sub(r3.Vec{X: t.n1.X, Y: t.n1.Y, Z: t.n1.Z},
		r3.Vec{X: t.n0.X, Y: t.n0.Y, Z: t.n0.Z})
	cos := r3.Cos(line, r3.Vec{X: 1.0, Y: 0.0, Z: 0.0})
	sin := r3.Cos(line, r3.Vec{X: 0.0, Y: 0.0, Z: 1.0})

	// Given the local xyz coordinate system	(x right, z downwards) and the global one (X right, Z
	// upwards), we can think of the correct coordinate transformation as two steps. First, the
	// rotation of this element around the local y-axis. And second, a rotation around the global X
	// axis to have the third axis point downwards. With s and c being sine and cosine of the angle
	// between element and local x-axis, both transformations look like this:
	//
	//                  1  0  0                   c  0  s                c   0   s
	// Around X-axis:   0 -1  0   Around y-axis:  0  1  0    Combined:   0  -1   0
	//                  0  0 -1                  -s  0  c                s   0  -c
	//
	// Note that y isn't relevant in a 2d element, but is listed for completeness. The y dimension is
	// left out for the rotation matrix below. Note that the resulting matrix is symmetric, but that
	// doesn't hold for all rotation matrices, so we still pick a non-symmetric matrix type.
	r := mat.NewDense(4, 4, nil)

	r.Set(0, 0, cos)
	r.Set(0, 1, sin)
	r.Set(1, 0, sin)
	r.Set(1, 1, -cos)

	r.Set(2, 2, cos)
	r.Set(2, 3, sin)
	r.Set(3, 2, sin)
	r.Set(3, 3, -cos)

	return r
}

func (t *truss2d) Indices(set map[xyz.Index]struct{}) {
	for _, index := range t.indicesAsArray() {
		set[index] = struct{}{}
	}
}

func (t *truss2d) NumNodes() uint {
	return 2
}

func (t *truss2d) ID() string {
	return t.id
}

// NewTruss2d returns a new linear-elastic 2d truss implementation, parametrised by its
// connectivity, cross section and Young's modulus (both assumed to be constant over the length),
// and hinges.
func NewTruss2d(
	id string,
	n0, n1 *xyz.Node,
	material *Material,
	disjoint Truss2dDisjoint) (Element, error) {
	// Any non-zero length is technically okay (even if not advisable), so use equality here to
	// detect actual zero-length trusses:
	if length(n0, n1) == 0.0 {
		return nil, fmt.Errorf("zero-length truss with nodes %v and %v", n0, n1)
	} else if id == "" {
		return nil, fmt.Errorf("element IDs can't be empty")
	}

	return &truss2d{id: id, n0: n0, n1: n1, material: material, disjoint: disjoint}, nil
}
