package deflect

import (
	"gonum.org/v1/gonum/spatial/r3"
)

func length(n0, n1 *Node) float64 {
	line := r3.Sub(r3.Vec{X: n1.X, Y: n1.Y, Z: n1.Z}, r3.Vec{X: n0.X, Y: n0.Y, Z: n0.Z})
	return r3.Norm(line)
}

func sineCosine2d(n0, n1 *Node) (sin, cos float64) {
	// Gonum's spatial/r2 package doesn't seem to have a way to create a matrix representation of a 2d
	// rotation (only for 3d), so we do it manually.
	line := r3.Sub(r3.Vec{X: n1.X, Y: n1.Y, Z: n1.Z}, r3.Vec{X: n0.X, Y: n0.Y, Z: n0.Z})
	cos = r3.Cos(line, r3.Vec{X: 1, Y: 0, Z: 0})
	sin = r3.Cos(line, r3.Vec{X: 0, Y: 0, Z: 1})

	return sin, cos
}

func directionCosine3d(n0, n1 *Node) (cxx, cxy, cxz float64) {
	line := r3.Sub(r3.Vec{X: n1.X, Y: n1.Y, Z: n1.Z}, r3.Vec{X: n0.X, Y: n0.Y, Z: n0.Z})
	cxx = r3.Cos(line, r3.Vec{X: 1, Y: 0, Z: 0})
	cxy = r3.Cos(line, r3.Vec{X: 0, Y: 1, Z: 0})
	cxz = r3.Cos(line, r3.Vec{X: 0, Y: 0, Z: 1})

	return cxx, cxy, cxz
}
