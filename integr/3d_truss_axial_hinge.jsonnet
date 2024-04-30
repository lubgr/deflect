local lib = import 'common.libsonnet';

local default = lib.Defaults();

local horizontal_axial_hinge(hinge) = {
  name: 'horizontal_axial_hinge',

  local l = 2.0,
  local F = 1.5e3,

  nodes: {
    A: [0, 5, -3],
    B: [-l, 5, -3],
    C: [-2 * l, 5, -3],
  },

  elements: {
    assert hinge == 'left' || hinge == 'right',
    AB: default.Truss3d(hinges=if hinge == 'right' then { B: ['Ux'] } else {}),
    BC: default.Truss3d(hinges=if hinge == 'left' then { B: ['Ux'] } else {}),
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: lib.Ux() + lib.Uy() + lib.Uz(),
    B: lib.Uz() + lib.Uy(),
    C: lib.Ux() + lib.Uy() + lib.Uz(),
  },

  neumann: {
    AB: lib.Fx(2 * F, l / 3),
    BC: lib.Fx(-F, 3 / 5 * l),
  },

  expected: {
    reaction: {
      A: lib.Fx(2 * F),
      C: lib.Fx(-F),
    },

    local EA = 30000e6 * 0.01,
    local uB = 2 * F / EA * 1 / 3 * l,
    local uC = F / EA * 2 / 5 * l,

    primary: {
      B: lib.Ux(if hinge == 'left' then -uB else uC),
    },

    interpolation: {
      AB: lib.Constant('Nx', 2 * F, range=[0, l / 3]) +
          lib.Constant('Nx', 0, range=[l / 3, l]) +
          lib.Linear('Ux', 0, uB, range=[0, l / 3]) +
          lib.Constant('Ux', uB, range=[l / 3, l]),
      BC: lib.Constant('Nx', 0, range=[0, 3 / 5 * l]) +
          lib.Constant('Nx', F, range=[3 / 5 * l, l]) +
          lib.Constant('Ux', -uC, range=[0, 3 / 5 * l]) +
          lib.Linear('Ux', -uC, 0, range=[3 / 5 * l, l]),
    },
  },
};

[
  horizontal_axial_hinge(hinge='left'),
  horizontal_axial_hinge(hinge='right'),
]
