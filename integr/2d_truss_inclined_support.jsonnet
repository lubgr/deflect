local lib = import 'common.libsonnet';

local default = lib.Defaults();

local triangle(angle) = {
  name: 'triangle_%d' % angle,
  description: 'Truss triangle with inclined support',

  local l = 1.0,
  local alpha = angle * lib.pi / 180.0,

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [0, 0, l],
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: default.Truss2d(),
    AC: default.Truss2d(),
    CB: default.Truss2d(),
  },

  dirichlet: {
    C: lib.Uz(),
    B: lib.Uz(),
  },

  local F = 10e3,

  neumann: {
    B: lib.Fx(F),
  },

  links: {
    A: lib.InclinedSupportUxUz(alpha),
  },

  expected: {
    reaction: {
      A: lib.Fx(-F) + lib.Fz(F / std.tan(alpha)),
      B: lib.Fz(0),
      C: lib.Fz(-F / std.tan(alpha)),
    },
  },
};

[
  triangle(20),
  triangle(45),
  triangle(90),  // Should be modelled with an ordinary Dirichlet BC in practice, but also works like this
  triangle(123.45678),
  triangle(315),
]
