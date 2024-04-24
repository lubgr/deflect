local lib = import 'common.libsonnet';

local default = lib.Defaults();

local triangle(H, V, a) = {
  name: 'triangle_%.2f_%.2f_%.2f' % [H, V, a],
  description: 'Truss triangle with H = %.2f, V = %.2f, a = %.2f' % [H, V, a],

  nodes: {
    A: [0, 0, 0],
    B: [2 * a, 0, 0],
    C: [a, 0, a],
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: default.Truss2d(),
    AC: default.Truss2d(),
    CB: default.Truss2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Uz(),
  },

  neumann: {
    C: lib.Fx(H) + lib.Fz(V),
  },

  expected: {
    reaction: {
      A: lib.Fx(-H) + lib.Fz((H + V) / (-2)),
      B: lib.Fz((H - V) / 2),
    },
    interpolation: {
      CB: lib.Constant('Nx', -(H - V) / std.sqrt(2)),
      AB: lib.Constant('Nx', (H - V) / 2),
      AC: lib.Constant('Nx', (H + V) / std.sqrt(2)),
    },
  },
};

[
  triangle(H=h, V=v, a=a)
  for h in [0, 10, -10]
  for v in [0, 10, -10]
  for a in [1, 5]
]
