package deflect

import (
	"fmt"
	"slices"
)

// oneDimElement holds state that is common to all 1d elements with 2 nodes; trusses, beams, and
// frames.
type oneDimElement struct {
	id       string
	n0, n1   *Node
	material *Material
	loads    []NeumannElementBC
}

func (e *oneDimElement) AddLoad(bc NeumannElementBC) bool {
	e.loads = append(e.loads, bc)
	return true
}

func (e *oneDimElement) RemoveLoad(bc NeumannElementBC) {
	for {
		i := slices.Index(e.loads, bc)

		if i == -1 {
			break
		}

		var zeroForGC NeumannElementBC
		e.loads[i] = zeroForGC
		e.loads = slices.Delete(e.loads, i, i+1)
	}
}

func (e *oneDimElement) NumNodes() uint {
	return 2
}

func (e *oneDimElement) ID() string {
	return e.id
}

func newOneDimElement(id string, n0, n1 *Node, material *Material) (oneDimElement, error) {
	// Any non-zero length is technically okay (even if not advisable), so use equality here to
	// detect actual zero-length trusses:
	if length(n0, n1) == 0.0 {
		return oneDimElement{}, fmt.Errorf("zero-length element with nodes %v and %v", n0, n1)
	} else if id == "" {
		return oneDimElement{}, fmt.Errorf("element IDs can't be empty")
	}

	return oneDimElement{id: id, n0: n0, n1: n1, material: material}, nil
}
