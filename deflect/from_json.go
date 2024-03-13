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
	// Disjoint maps node IDs to a sequence of string representations of degrees of freedoms that are
	// hinged. Example: {"A": ["Ux"], "B", ["Uz", "Phiy"]}.
	Disjoint map[string][]string
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

	nodes := translateNodes(tmp.Nodes)
	lookup := matLookupFct(materials, cs)
	elements, errElements := translateElements(tmp.Elements, nodes, lookup)
	if errElements != nil {
		return Problem{}, fmt.Errorf(
			"couldn't construct element instances: %w",
			errElements,
		)
	}

	dirichletBCs, errDirichlet := translateDirichletBCs(tmp.Dirichlet, nodes)
	neumannNodalBCs, errNeumann := translateNodalNeumannBCs(tmp.Neumann, nodes)
	links, linkBCs, errLinks := translateAngularLinks(tmp.Links, nodes)

	if err := errors.Join(errDirichlet, errNeumann, errLinks); err != nil {
		return Problem{}, fmt.Errorf("couldn't construct BCs: %w", err)
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

func translateNodes(from map[string][3]float64) []Node {
	nodes := make([]Node, 0, len(from))
	for id, coor := range from {
		nodes = append(nodes, Node{ID: id, X: coor[0], Y: coor[1], Z: coor[2]})
	}
	return nodes
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

	for id, desc := range from {
		mat, errLookup := material(desc.Material, desc.CS)
		if errLookup != nil {
			return elements, errLookup
		}

		if desc.Kind == "truss2d" {
			n0, n1, errNodes := determineNodes(id, desc.Nodes, nodes)
			disjoint, errDisjoint := determineDisjoints(desc.Disjoint, n0.ID, n1.ID)
			truss, errCreate := NewTruss2d(id, n0, n1, &mat, disjoint)

			if err := errors.Join(errNodes, errDisjoint, errCreate); err != nil {
				return elements, err
			}

			elements = append(elements, truss)
		}
	}

	return elements, nil
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

func determineDisjoints(from map[string][]string, n0, n1 string) (Truss2dDisjoint, error) {
	idx0 := slices.Index(from[n0], "Ux")
	idx1 := slices.Index(from[n1], "Ux")

	switch {
	case idx0 != -1 && idx1 != -1:
		return Truss2dDisjointNone,
			fmt.Errorf("only a single Ux entry makes sense in a 2d truss disjoint description")
	case idx0 != -1:
		return Truss2dDisjointFirstNode, nil
	case idx1 != -1:
		return Truss2dDisjointSecondNode, nil
	default:
		return Truss2dDisjointNone, nil
	}
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
	dofs := map[string]Dof{
		"Ux":   Ux,
		"Uz":   Uz,
		"Uy":   Uy,
		"Phiy": Phiy,
		"Phiz": Phiz,
		"Phix": Phix,
	}
	result := make([]NodalValue, 0, len(from))
	var err error

	for nodeID, dofDescriptions := range from {
		for _, dofToValue := range dofDescriptions {
			for dofName, value := range dofToValue {
				dof, ok := dofs[dofName]

				if !ok {
					err = errors.Join(err, fmt.Errorf("unknown degree of freedom '%v'", dofName))
					continue
				}

				index := Index{NodalID: nodeID, Dof: dof}
				result = append(result, NodalValue{Index: index, Value: value})
			}
		}
	}

	return result, err
}

func translateNodalNeumannBCs(
	from map[string][]neumannDescription,
	nodes []Node,
) ([]NodalValue, error) {
	dofs := map[string]Dof{
		"Fx": Ux,
		"Fz": Uz,
		"Fy": Uy,
		"My": Phiy,
		"Mz": Phiz,
		"Mx": Phix,
	}
	result := make([]NodalValue, 0, len(from))
	var err error

	for nodeID, neumannDesc := range from {
		for _, desc := range neumannDesc {
			// When the object describes an element Neumann BC, the sequence is empty.
			for name, value := range desc.Nodal {
				dof, ok := dofs[name]

				if !ok {
					err = errors.Join(err, fmt.Errorf("unknown Neumann BC type '%v'", name))
					continue
				}

				index := Index{NodalID: nodeID, Dof: dof}
				result = append(result, NodalValue{Index: index, Value: value})
			}
		}
	}

	return result, err
}

func translateAngularLinks(
	from map[string][]dirichletAngularLink,
	nodes []Node,
) (links []Transformer, bcs []NodalValue, err error) {
	dofs := map[string]Dof{
		"Ux": Ux,
		"Uz": Uz,
		"Uy": Uy,
	}
	links = make([]Transformer, 0, len(from))
	bcs = make([]NodalValue, 0, len(from))

	for nodeID, linkDesc := range from {
		for _, desc := range linkDesc {
			dofFrom, okFrom := dofs[desc.From]
			dofTo, okTo := dofs[desc.To]

			if !(okFrom && okTo) {
				err = errors.Join(
					err,
					fmt.Errorf("unknown angular displacement BC type '%v' or '%v'", desc.From, desc.To),
				)
				continue
			}

			fromIndex := Index{NodalID: nodeID, Dof: dofFrom}
			toIndex := Index{NodalID: nodeID, Dof: dofTo}
			transformer, bc := NewInclinedSupport(fromIndex, toIndex, desc.Angle)
			links = append(links, transformer)
			bcs = append(bcs, bc)
		}
	}

	return links, bcs, err
}
