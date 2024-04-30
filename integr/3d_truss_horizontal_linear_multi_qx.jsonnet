local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

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

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: bvp.Uz() + bvp.Uy(),
    B: bvp.Ux() + bvp.Uy() + bvp.Uz(),
  },

  neumann: {
    AB: bvp.qx([0, 2e3]) + bvp.qx([-2e3, 5e3]) + bvp.qx([-2e3, 0]) + bvp.qx([2e3, -5e3]),
  },

  expected: {
    reaction: {
      B: test.Fx(0),
    },
    interpolation: {
      local Nx(x) = 2e3 * x - 2e3 / l * x * x,
      AB: test.Quadratic('Nx', eval=test.Samples(Nx, 0, l, 7)),
    },
  },
}
