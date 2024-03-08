package deflect

import "fmt"

// oneDimElement holds state that is common to all 1d elements with 2 nodes; trusses, beams, and
// frames.
type oneDimElement struct {
	id       string
	n0, n1   *Node
	material *Material
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
