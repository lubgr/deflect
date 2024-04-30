local lib = import 'common.libsonnet';

local default = lib.Defaults();

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

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: lib.Uz() + lib.Uy(),
    B: lib.Uz() + lib.Uy(),
    C: lib.Ux() + lib.Uy() + lib.Uz(),
  },

  neumann: {
    A: lib.Fx(F),
  },

  expected: {
    local uA = F / (30000e6 * 0.01) * 2 * l,
    local uB = uA / 2,

    reaction: {
      C: lib.Fx(-F),
    },
    primary: {
      A: lib.Ux(uA),
      B: lib.Ux(uB),
    },
    interpolation: {
      // The Ux assertions ensure that the right sign is used for the nodal Ux value
      // when integrating ε to become the displacement field.
      AB: lib.Constant('Nx', -F) +
          lib.Linear('Ux', uA, uB),
      BC: lib.Constant('Nx', -F) +
          lib.Linear('Ux', uB, 0),
    },
  },
}
