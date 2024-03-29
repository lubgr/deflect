# Choose Global and Local CO-System for Mesh and Elements

## Context and Problem Statement

The solver base needs two central coordinate systems: the local coordinate system for element
formulations and the global coordinate system. The global one will be present everywhere, e.g. in
tests or when element formulations transform quantities from their local coordinate system into the
global one. There is no way to parametrise this, since it would have to be parametrised in terms of
another root coordinate system. The local coordinate system will be present in element formulations
and related functionality like cross section characteristics. For obvious reasons, we stick with the
same local coordinate system for all classical 1d structural elements like beams and trusses.

This choice is independent from the coordinate system of a UI. It would be the simplest approach to
mirror these coordinate systems in a UI, but having transformations in between is still an option.

## Considered Options

All right-handed coordinate system can be considered for both local and global coordinate systems.
It could make sense to try and determine what has most consensus across text books (but that's not a
fun task, so we don't).

## Decision Outcome

The choice is somewhat arbitrary and according to my intuition at the time of writing.

The global x-axis is horizontal, pointing from left to right. The global z-axis points upwards. The
orientation of the global y-axis follows from the right hand rule; if the basis vectors were drawn
on a sheet of paper, the global y-axis would point "into" the paper.

The element-local x-axis is following the one dimension of classical 1d elements like beams, frames,
and trusses. It points from the first node to the second node. For an x-axis from left to right
(e.g., a two-node element with the first node on the left and the second node on the right), the
local z-axis points downwards. The orientation of the local y-axis again follows from that; on the
hypothetical piece of paper, the local y-axis points "out of" the paper in this case.

This setup can be different to what is found in textbooks. Formulae for e.g. linear-elastic beam
stiffness matrices when transformed to global coordinates are then different from what we implement
here, specifically the signs of the transformation matrices.

To illustrate a 2d coordinate transformation from local element coordinates to global coordinates,
the transformation can be thought of as two steps. First, the rotation of an element around the
local y-axis. And second, a rotation around the global X axis to let the third axis point downwards.
With s and c being sine and cosine of the angle between element-local x-axis and global X-axis, both
transformations look like this:

                     1  0  0                   c  0  s                c   0   s
    Around X-axis:   0 -1  0   Around y-axis:  0  1  0    Combined:   0  -1   0
                     0  0 -1                  -s  0  c                s   0  -c
