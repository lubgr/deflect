package deflect

import "testing"

var exampleMat = Material{
	CrossSection: &rectangular{b: 0.1, h: 0.1, roll: 0},
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

func TestOneDimElementBCStorage(t *testing.T) {
	n0, n1 := oneDimElementTestNodes()
	e, _ := newOneDimElement("ID", n0, n1, &exampleMat)
	data := []int{100, 200, 300, 400, 500}
	load := []NeumannElementBC{&data[0], &data[1], &data[2], &data[3], &data[4]}

	e.AddLoad(load[0])
	e.AddLoad(load[1])
	e.AddLoad(load[2])
	e.AddLoad(load[2]) // Adding the same load multiple times is supported
	e.AddLoad(load[3])
	e.AddLoad(load[2])

	if n := len(e.loads); n != 6 {
		t.Errorf("Each Element BC must be accounted for once, expected 6 BC, got %v", n)
	}

	removal := []struct {
		remove    NeumannElementBC
		remaining int
	}{
		{remove: load[4], remaining: 6},
		{remove: load[0], remaining: 5},
		{remove: load[3], remaining: 4},
		{remove: load[2], remaining: 1},
		{remove: load[1], remaining: 0},
		{remove: load[1], remaining: 0},
	}

	for _, test := range removal {
		e.RemoveLoad(test.remove)
		if n := len(e.loads); n != test.remaining {
			t.Errorf("Expected %v remaining BC after preceding removal, got %v", test.remaining, n)
		}
	}
}
