# Use Symbolic Indices to Global Matrix/Vector Indices

## Context and Problem Statement

Local element matrices and vectors have to be connected to global ones by establishing a mapping
between local and global matrix/vector indices. The same applies to nodal Dirichlet and Neumann
boundary conditions, e.g., what is the matrix index for a Dirichlet boundary condition that
prescribes the vertical displacement at node "A"? How and where is this mapping implemented?

## Considered Options

- "ID array" as covered in textbooks
- Custom hash map with symbolic indices

## Decision Outcome

We go with a hash map of symbolic indices to integral matrix indices. A "symbolic" or "semantic"
index is composed of a nodal ID and a degree of freedom identifier. `{"A", Ux}` refers to the
horizontal displacement at node "A" and in global coordinates. Element quantities and boundary
conditions are constructed using such semantic indices, i.e., the carry these symbolic indices as
part of their object state.

It is the responsibility of each element formulation to write its contribution directly to the
_global_ tangent and residual. A lookup table is provided to this end. This table is created as one
of the first steps of solving a boundary value problem, and it maps symbolic indices to integral
indices into the global tangent matrix and residual vector. For example, a 2d truss element that is
connected to nodes "C" and "E", constructs four such symbolic indices: `{"C", Ux}`, `{"C", Uz}`,
`{"E", Ux}`, and `{"E", Uz}`. Using the mapping, these could translate into plain matrix indices 5,
6, 90, 91.

Adding values directly to the global tangent and residual is more efficient than constructing
intermediate local matrices (although an element might does that as an implementation detail). It
does give elements more responsibility, and they can theoretically mess up global matrices instead
of just local ones. However, such bugs would invalidate the system of equations to be solved anyhow.

The mapping can also be used to partition the global system of equations into a part with
Dirichlet-constrained and free degrees of freedom. This allows us to solve for the free degrees of
freedom by using a subset of the global matrices without manually rearranging rows and columns or
copying coefficients into different matrices.

This solution naturally allows for heterogeneous degrees of freedom per node. For example, in a
boundary value problem with both truss and frame elements, there can be nodes without rotational
degrees of freedom (when no frame connects to this node); no storage is wasted for representing the
unused rotational degree of freedom at these nodes.
