local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local triangle(angle) = {
  name: 'triangle_%d' % angle,
  description: 'Truss triangle with inclined support',

  local l = 1.0,
  local alpha = angle * bvp.pi / 180.0,

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [0, 0, l],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: bvp.Truss2d(),
    AC: bvp.Truss2d(),
    CB: bvp.Truss2d(),
  },

  dirichlet: {
    C: bvp.Uz(),
    B: bvp.Uz(),
  },

  local F = 10e3,

  neumann: {
    B: bvp.Fx(F),
  },

  links: {
    A: bvp.InclinedSupportUxUz(alpha),
  },

  expected: {
    reaction: {
      A: test.Fx(-F) + test.Fz(F / std.tan(alpha)),
      B: test.Fz(0),
      C: test.Fz(-F / std.tan(alpha)),
    },
    interpolation: {
      AB: test.Constant('Nx', F),
      CB: test.Constant('Nx', 0),
      AC: test.Constant('Nx', -F / std.tan(alpha)),
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
