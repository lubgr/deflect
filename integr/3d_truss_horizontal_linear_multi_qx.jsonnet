local lib = import 'common.libsonnet';

local default = lib.Defaults();

{
  name: 'horizontal_linear_multi_nx',

  local l = 3.0,

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  elements: {
    AB: default.Truss3d(),
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: lib.Uz() + lib.Uy(),
    B: lib.Ux() + lib.Uy() + lib.Uz(),
  },

  neumann: {
    AB: lib.qx([0, 2e3]) + lib.qx([-2e3, 5e3]) + lib.qx([-2e3, 0]) + lib.qx([2e3, -5e3]),
  },

  expected: {
    reaction: {
      B: lib.Fx(0),
    },
    interpolation: {
      local Nx(x) = 2e3 * x - 2e3 / l * x * x,
      AB: lib.Quadratic('Nx', eval=lib.Samples(Nx, 0, l, 7)),
    },
  },
}
