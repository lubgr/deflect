local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

{
  name: 'horizontal_two_elements',

  local l = 2.0,
  local F = 1.5e3,

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [2 * l, 0, 0],
  },

  elements: {
    AB: default.Truss3d(),
    BC: default.Truss3d(),
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: bvp.Uz() + bvp.Uy(),
    B: bvp.Uz() + bvp.Uy(),
    C: bvp.Ux() + bvp.Uy() + bvp.Uz(),
  },

  neumann: {
    A: bvp.Fx(F),
  },

  expected: {
    local uA = F / (30000e6 * 0.01) * 2 * l,
    local uB = uA / 2,

    reaction: {
      C: test.Fx(-F),
    },
    primary: {
      A: test.Ux(uA),
      B: test.Ux(uB),
    },
    interpolation: {
      // The Ux assertions ensure that the right sign is used for the nodal Ux value
      // when integrating ε to become the displacement field.
      AB: test.Constant('Nx', -F) +
          test.Linear('Ux', uA, uB),
      BC: test.Constant('Nx', -F) +
          test.Linear('Ux', uB, 0),
    },
  },
}
