package elmt

import (
	"github.com/lubgr/deflect/xyz"
	"gonum.org/v1/gonum/spatial/r3"
)

func length(n0, n1 *xyz.Node) float64 {
	line := r3.Sub(r3.Vec{X: n1.X, Y: n1.Y, Z: n1.Z}, r3.Vec{X: n0.X, Y: n0.Y, Z: n0.Z})
	return r3.Norm(line)
}
