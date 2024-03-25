package deflect

import (
	"encoding/json"
	"errors"
	"fmt"
	"slices"
)

type matDescription struct {
	Kind      string
	Parameter map[string]float64
}

type csDescription struct {
	Kind      string
	Parameter map[string]float64
}

type nodalValues map[string]float64

type dirichletAngularLink struct {
	From  string
	To    string
	Angle float64
}

// neumannDescription is a simplistic sum type, where nothing enforces that is isn't used as a
// product type. Should be ok, as it's used in a very small scope.
type neumannDescription struct {
	// When Nodal is populated, a nodal Neumann BC is described.
	Nodal map[string]float64
	// When Nodal is not populated, we have an element Neumann BC.
	Element struct {
		Kind     string
		Degree   string
		Values   []float64
		Position []float64
	}
}

func (nd *neumannDescription) UnmarshalJSON(data []byte) error {
	if string(data) == "null" || len(data) == 0 {
		return nil
	}

	var tmp map[string]any
	if err := json.Unmarshal(data, &tmp); err != nil {
		return fmt.Errorf("failed to unmarshal Neumann BC JSON: %w", err)
	}

	_, elementHint1 := tmp["kind"].(string)
	_, elementHint2 := tmp["degree"].(string)

	if elementHint1 && elementHint2 {
		return json.Unmarshal(data, &nd.Element)
	}

	return json.Unmarshal(data, &nd.Nodal)
}

type elmtDescription struct {
	Kind     string
	Nodes    []string
	CS       string
	Material string
	// Hinge maps node IDs to a sequence of string representations of degrees of freedoms that are
	// hinged. Example: {"A": ["Ux"], "B", ["Uz", "Phiy"]}.
	Hinges map[string][]string
}

// FromJSON parses the given JSON data and constructs a boundary value problem from it.
func FromJSON(data []byte) (Problem, error) {
	var tmp struct {
		Nodes        map[string][3]float64
		Material     map[string]matDescription
		Crosssection map[string]csDescription
		Elements     map[string]elmtDescription
		Dirichlet    map[string][]nodalValues
		Links        map[string][]dirichletAngularLink
		Neumann      map[string][]neumannDescription
	}

	if err := json.Unmarshal(data, &tmp); err != nil {
		return Problem{}, fmt.Errorf("couldn't build Bvp from JSON: %w", err)
	}

	cs, errCS := translateCrossSections(tmp.Crosssection)
	materials, errMaterials := translateMaterials(tmp.Material)

	if err := errors.Join(errCS, errMaterials); err != nil {
		return Problem{}, fmt.Errorf("couldn't build CS/material maps: %w", err)
	}

	nodes, errNodes := translateNodes(tmp.Nodes)
	if errNodes != nil {
		return Problem{}, fmt.Errorf("construct mesh: %w", errNodes)
	}

	elements, errElements := translateElements(tmp.Elements, nodes, matLookupFct(materials, cs))
	if errElements != nil {
		return Problem{}, fmt.Errorf("construct mesh: %w", errElements)
	}

	dirichletBCs, errDirichlet := translateDirichletBCs(tmp.Dirichlet, nodes)
	neumannNodalBCs, errNeumann0 := translateNodalNeumannBCs(tmp.Neumann, nodes)
	errNeumann1 := translateAndApplyElementNeumannBCs(tmp.Neumann, elements)
	links, linkBCs, errLinks := translateAngularLinks(tmp.Links, nodes)

	if err := errors.Join(errDirichlet, errNeumann0, errNeumann1, errLinks); err != nil {
		return Problem{}, fmt.Errorf("construct BCs: %w", err)
	} else if len(dirichletBCs)+len(linkBCs) == 0 {
		return Problem{}, errors.New("can't construct a BVP with no Dirichlet BC")
	}

	result := Problem{
		Nodes:        nodes,
		Elements:     elements,
		Dirichlet:    append(dirichletBCs, linkBCs...),
		Neumann:      neumannNodalBCs,
		EqTransforms: links,
	}

	return result, nil
}

func translateNodes(from map[string][3]float64) ([]Node, error) {
	nodes := make([]Node, 0, len(from))
	for id, coor := range from {
		nodes = append(nodes, Node{ID: id, X: coor[0], Y: coor[1], Z: coor[2]})
	}

	if len(nodes) == 0 {
		return nodes, errors.New("no nodes created")
	}

	return nodes, nil
}

func matLookupFct(
	mats map[string]LinearElastic,
	css map[string]CrossSection,
) func(string, string) (Material, error) {
	return func(mat, cs string) (Material, error) {
		result := Material{}
		var okMat, okCS bool

		result.LinearElastic, okMat = mats[mat]
		result.CrossSection, okCS = css[cs]

		if !okMat {
			return result, fmt.Errorf("material '%v' not found", mat)
		}
		if !okCS {
			return result, fmt.Errorf("cross section '%v' not found", cs)
		}

		return result, nil
	}
}

func translateElements(
	from map[string]elmtDescription,
	nodes []Node,
	material func(string, string) (Material, error),
) ([]Element, error) {
	elements := make([]Element, 0, len(from))
	var err error

	for id, desc := range from {
		mat, errLookup := material(desc.Material, desc.CS)
		if errLookup != nil {
			err = errors.Join(fmt.Errorf("couldn't lookup material for element %v: %w", id, errLookup))
			continue
		}

		switch desc.Kind {
		case "truss2d":
			n0, n1, errNodes := determineNodes(id, desc.Nodes, nodes)
			hinges, errHinges := formHingeMap(desc.Hinges, n0.ID, n1.ID)
			truss, errCreate := NewTruss2d(id, n0, n1, &mat, hinges)

			if errJoined := errors.Join(errNodes, errHinges, errCreate); errJoined != nil {
				err = errors.Join(err, fmt.Errorf("create 2d truss '%v': %w", id, errJoined))
				continue
			}

			elements = append(elements, truss)
		case "frame2d":
			n0, n1, errNodes := determineNodes(id, desc.Nodes, nodes)
			hinges, errHinges := formHingeMap(desc.Hinges, n0.ID, n1.ID)
			frame, errCreate := NewFrame2d(id, n0, n1, &mat, hinges)

			if errJoined := errors.Join(errNodes, errHinges, errCreate); errJoined != nil {
				err = errors.Join(err, fmt.Errorf("create 2d frame '%v': %w", id, errJoined))
				continue
			}

			elements = append(elements, frame)
		}
	}

	if len(elements) == 0 {
		return elements, errors.New("no elements created")
	}

	return elements, err
}

// determineNodes looks up pointers to nodes from the given slice of nodes, using either nodeIDs
// when it's of length 2 or the elmtID string, which is then expected to be a two-character string.
// Valid examples:
// - elmtID is "DF", nodeIDs is [], and there are {"D", x1, y1, z1} and {"F", x2, y2, z2} in nodes
// - elmtID is "1234", nodeIDs is ["D", "F"], and there are {"D", ...} and {"F", ...} in nodes
// Note that we could be more efficient by maintaining the nodes in a hash map using their IDs. We
// should do this once this linear lookup shows up in a profile.
func determineNodes(
	elmtID string,
	nodeIDs []string,
	nodes []Node,
) (*Node, *Node, error) {
	if len(nodeIDs) == 2 {
		n0, err0 := scanForNode(nodeIDs[0], nodes)
		n1, err1 := scanForNode(nodeIDs[1], nodes)
		return n0, n1, errors.Join(err0, err1)
	} else if len(nodeIDs) > 0 {
		return nil, nil, fmt.Errorf("expected node ID slice of length 0 or 2, got %v IDs", len(nodeIDs))
	}

	if len(elmtID) != 2 {
		return nil, nil, fmt.Errorf(
			"can't infer connected node IDs from element ID '%v' (need 2 chars)",
			elmtID,
		)
	}

	n0, n1 := string(elmtID[0]), string(elmtID[1])
	return determineNodes(elmtID, []string{n0, n1}, nodes)
}

func scanForNode(nodeID string, nodes []Node) (*Node, error) {
	idx := slices.IndexFunc(nodes, func(n Node) bool { return n.ID == nodeID })
	if idx == -1 {
		return nil, fmt.Errorf("nodal ID '%v' not found in sequence of %v nodes", nodeID, len(nodes))
	}
	return &nodes[idx], nil
}

func formHingeMap(from map[string][]string, n0, n1 string) (map[Index]struct{}, error) {
	dofs := dofLookup{
		context: "construct hinge",
		dofs: map[string]Dof{
			"Ux":   Ux,
			"Uz":   Uz,
			"Uy":   Uy,
			"Phiy": Phiy,
			"Phiz": Phiz,
			"Phix": Phix,
		}}
	hinges := map[Index]struct{}{}
	var err error

	for nodeID, dofNames := range from {
		if nodeID != n0 && nodeID != n1 {
			err = errors.Join(
				err,
				fmt.Errorf("hinge node %v must be the element's nodes (%v or %v)", nodeID, n0, n1),
			)
			continue
		}

		for _, dofName := range dofNames {
			dof, ok := dofs.Lookup(dofName)
			if !ok {
				continue
			}
			hinges[Index{NodalID: nodeID, Dof: dof}] = struct{}{}
		}
	}

	return hinges, dofs.FinaliseJoin(err)
}

func translateMaterials(from map[string]matDescription) (map[string]LinearElastic, error) {
	materials := map[string]LinearElastic{}

	for id, desc := range from {
		E, found0 := desc.Parameter["E"]
		nu, found1 := desc.Parameter["nu"]
		rho, found2 := desc.Parameter["rho"]

		if !found0 || !found1 || !found2 {
			return materials, fmt.Errorf("material parameters 'E', 'nu', and/or 'rho' not found")
		}

		materials[id] = LinearElastic{
			YoungsModulus: E,
			PoissonsRatio: nu,
			Density:       rho,
		}
	}

	return materials, nil
}

func translateCrossSections(from map[string]csDescription) (map[string]CrossSection, error) {
	cs := map[string]CrossSection{}
	var errs []error
	fail := func(msg string, a ...any) { errs = append(errs, fmt.Errorf(msg, a...)) }

	for id, desc := range from {
		switch desc.Kind {

		case "rectangle":
			b, foundB := desc.Parameter["b"]
			h, foundH := desc.Parameter["h"]

			if !foundB || !foundH {
				fail("cross section parameters 'b', and/or 'h' not found")
				continue
			} else if rect, err := NewRectangularCrossSection(b, h); err != nil {
				fail("create rectangular cross section: %w", err)
			} else {
				cs[id] = rect
			}
		case "constants":
			constants, err := NewConstantsCrossSections(desc.Parameter)

			if err != nil {
				fail("create constants-based cross section: %w", err)
			} else {
				cs[id] = constants
			}
		default:
			fail("unsupported cross section '%v'", desc.Kind)
		}
	}

	return cs, errors.Join(errs...)
}

func translateDirichletBCs(
	from map[string][]nodalValues,
	nodes []Node,
) ([]NodalValue, error) {
	dofs := dofLookup{
		context: "construct Dirichlet BC",
		dofs: map[string]Dof{
			"Ux":   Ux,
			"Uz":   Uz,
			"Uy":   Uy,
			"Phiy": Phiy,
			"Phiz": Phiz,
			"Phix": Phix,
		}}
	result := make([]NodalValue, 0, len(from))

	for nodeID, dofDescriptions := range from {
		for _, dofToValue := range dofDescriptions {
			for dofName, value := range dofToValue {
				dof, ok := dofs.Lookup(dofName)
				if !ok {
					continue
				}

				index := Index{NodalID: nodeID, Dof: dof}
				result = append(result, NodalValue{Index: index, Value: value})
			}
		}
	}

	return result, dofs.FinaliseJoin(nil)
}

func translateNodalNeumannBCs(
	from map[string][]neumannDescription,
	nodes []Node,
) ([]NodalValue, error) {
	dofs := dofLookup{
		context: "construct nodal Neumann BC",
		dofs: map[string]Dof{
			"Fx": Ux,
			"Fz": Uz,
			"Fy": Uy,
			"My": Phiy,
			"Mz": Phiz,
			"Mx": Phix,
		}}
	result := make([]NodalValue, 0, len(from))

	for nodeID, neumannDesc := range from {
		for _, desc := range neumannDesc {
			// When the object describes an element Neumann BC, the sequence is empty.
			for dofName, value := range desc.Nodal {
				dof, ok := dofs.Lookup(dofName)
				if !ok {
					continue
				}

				index := Index{NodalID: nodeID, Dof: dof}
				result = append(result, NodalValue{Index: index, Value: value})
			}
		}
	}

	return result, dofs.FinaliseJoin(nil)
}

func translateAndApplyElementNeumannBCs(
	from map[string][]neumannDescription,
	elements []Element,
) error {
	var err error

	for elmtID, neumannDesc := range from {
		for _, desc := range neumannDesc {
			if len(desc.Nodal) != 0 {
				// This describes a nodal Neumann BC, which is handled elsewhere.
				continue
			}

			elmt, errElmt := scanForElement(elmtID, elements)
			if errElmt != nil {
				err = errors.Join(err, fmt.Errorf("find element for applying Neumann BC: %w", errElmt))
				continue
			}

			load, errLoad := translateElementNeumannBC(&desc)

			if errLoad != nil {
				err = errors.Join(err, fmt.Errorf("construct load for element %v: %w", elmtID, errLoad))
				continue
			} else if !elmt.AddLoad(load) {
				err = errors.Join(err, fmt.Errorf("couldn't apply load on element %v", elmtID))
			}
		}
	}

	return err
}

func translateElementNeumannBC(desc *neumannDescription) (NeumannElementBC, error) {
	dofs := dofLookup{
		context: "construct element Neumann BC",
		dofs: map[string]Dof{
			"Fx": Ux,
			"Fz": Uz,
			"Fy": Uy,
			"My": Phiy,
			"Mz": Phiz,
			"Mx": Phix,
			"qx": Ux,
			"qz": Uz,
			"qy": Uy,
			"my": Phiy,
			"mz": Phiz,
			"mx": Phix,
		}}
	var load NeumannElementBC

	kind, ok := dofs.Lookup(desc.Element.Kind)

	if !ok {
		return load, dofs.FinaliseJoin(nil)
	}

	values, pos := desc.Element.Values, desc.Element.Position
	numValues, numPos := len(values), len(pos)
	var err error

	switch desc.Element.Degree {
	case "dirac":
		if numValues != 1 || numPos != 1 {
			err = fmt.Errorf("concentrated load needs single value/position, not %v and %v", values, pos)
		} else {
			load, err = NewElementConcentratedLoad(kind, pos[0], values[0])
		}
	case "constant":
		if numValues != 1 {
			err = fmt.Errorf("constant load must have 1 value, got %v", numValues)
		} else {
			load = NewElementConstantLoad(kind, values[0])
		}
	case "linear":
		if numValues != 2 {
			err = fmt.Errorf("linear load must have 2 values, got %v", numValues)
		} else {
			load = NewElementLinearLoad(kind, values[0], values[1])
		}
	default:
		err = fmt.Errorf("unknown polynomial degree for Element BC '%v'", desc.Element.Degree)
	}

	return load, err
}

func scanForElement(elmtID string, elements []Element) (Element, error) {
	idx := slices.IndexFunc(elements, func(e Element) bool { return e.ID() == elmtID })
	if idx == -1 {
		return nil, fmt.Errorf(
			"element ID '%v' not found in sequence of %v elements",
			elmtID,
			len(elements),
		)
	}
	return elements[idx], nil
}

func translateAngularLinks(
	from map[string][]dirichletAngularLink,
	nodes []Node,
) (links []Transformer, bcs []NodalValue, err error) {
	dofs := dofLookup{
		context: "construct inclined support",
		dofs: map[string]Dof{
			"Ux": Ux,
			"Uz": Uz,
			"Uy": Uy,
		}}
	links = make([]Transformer, 0, len(from))
	bcs = make([]NodalValue, 0, len(from))

	for nodeID, linkDesc := range from {
		for _, desc := range linkDesc {
			dofFrom, okFrom := dofs.Lookup(desc.From)
			dofTo, okTo := dofs.Lookup(desc.To)

			if !(okFrom && okTo) {
				continue
			}

			fromIndex := Index{NodalID: nodeID, Dof: dofFrom}
			toIndex := Index{NodalID: nodeID, Dof: dofTo}
			transformer, bc := NewInclinedSupport(fromIndex, toIndex, desc.Angle)
			links = append(links, transformer)
			bcs = append(bcs, bc)
		}
	}

	return links, bcs, dofs.FinaliseJoin(nil)
}

type dofLookup struct {
	dofs    map[string]Dof
	context string
	err     error
}

func (dl *dofLookup) Lookup(key string) (Dof, bool) {
	dof, ok := dl.dofs[key]

	if !ok {
		dl.err = errors.Join(
			dl.err,
			fmt.Errorf("%v: '%v' not mapped to degree of freedom", dl.context, key),
		)
	}

	return dof, ok
}

func (dl *dofLookup) FinaliseJoin(other error) error {
	return errors.Join(other, dl.err)
}
