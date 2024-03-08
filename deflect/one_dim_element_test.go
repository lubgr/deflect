package deflect

import "testing"

var exampleMat = Material{
	CrossSection: &rectangular{0.1, 0.1},
	LinearElastic: LinearElastic{
		YoungsModulus: 100000.0 / (1e-3 * 1e-3),
		PoissonsRatio: 0.3,
		Density:       1,
	},
}

func oneDimElementTestNodes() (n0, n1 *Node) {
	n0 = &Node{ID: "A", X: 0, Z: 0}
	n1 = &Node{ID: "B", X: 1, Z: 0}
	return
}

func TestOneDimElementCtorFailsOnZeroLength(t *testing.T) {
	n0 := &Node{ID: "A", X: 2, Y: 3, Z: -5}
	n1 := &Node{ID: "B", X: 2, Y: 3, Z: -5}

	actual, err := newOneDimElement("AB", n0, n1, &exampleMat)

	if err == nil {
		t.Errorf("1d element base creation succeeded despite zero length, got %v", actual)
	}
}

func TestOneDimElementCtorFailsOnEmptyElementId(t *testing.T) {
	n0, n1 := oneDimElementTestNodes()
	actual, err := newOneDimElement("", n0, n1, &exampleMat)

	if err == nil {
		t.Errorf("1d element creation succeeded despite empty element id, got %v", actual)
	}
}

func TestOneDimElementReturnsID(t *testing.T) {
	// Stupid test, candidate for removal from day one...
	n0, n1 := oneDimElementTestNodes()
	id := "ID123"
	element, _ := newOneDimElement(id, n0, n1, &exampleMat)

	if actual := element.ID(); actual != id {
		t.Errorf("Expected element to return an ID as %v, but got %v", id, actual)
	}
}

func TestOneDimElementNumNodes(t *testing.T) {
	n0, n1 := oneDimElementTestNodes()
	element, _ := newOneDimElement("ID", n0, n1, &exampleMat)
	actual := element.NumNodes()

	if actual != 2 {
		t.Errorf("Expected 1d element to indicate 2 connected nodes, got %v", actual)
	}
}
