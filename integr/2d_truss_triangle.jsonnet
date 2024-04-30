local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local triangle(H, V, a) = {
  name: 'triangle_%.2f_%.2f_%.2f' % [H, V, a],
  description: 'Truss triangle with H = %.2f, V = %.2f, a = %.2f' % [H, V, a],

  nodes: {
    A: [0, 0, 0],
    B: [2 * a, 0, 0],
    C: [a, 0, a],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: bvp.Truss2d(),
    AC: bvp.Truss2d(),
    CB: bvp.Truss2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    B: bvp.Uz(),
  },

  neumann: {
    C: bvp.Fx(H) + bvp.Fz(V),
  },

  expected: {
    reaction: {
      A: test.Fx(-H) + test.Fz((H + V) / (-2)),
      B: test.Fz((H - V) / 2),
    },
    interpolation: {
      CB: test.Constant('Nx', -(H - V) / std.sqrt(2)),
      AB: test.Constant('Nx', (H - V) / 2),
      AC: test.Constant('Nx', (H + V) / std.sqrt(2)),
    },
  },
};

[
  triangle(H=h, V=v, a=a)
  for h in [0, 10, -10]
  for v in [0, 10, -10]
  for a in [1, 5]
]
