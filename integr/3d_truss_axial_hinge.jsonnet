local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

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
    AB: bvp.Truss3d(hinges=if hinge == 'right' then { B: ['Ux'] } else {}),
    BC: bvp.Truss3d(hinges=if hinge == 'left' then { B: ['Ux'] } else {}),
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: bvp.Ux() + bvp.Uy() + bvp.Uz(),
    B: bvp.Uz() + bvp.Uy(),
    C: bvp.Ux() + bvp.Uy() + bvp.Uz(),
  },

  neumann: {
    AB: bvp.Fx(2 * F, l / 3),
    BC: bvp.Fx(-F, 3 / 5 * l),
  },

  expected: {
    reaction: {
      A: test.Fx(2 * F),
      C: test.Fx(-F),
    },

    local EA = 30000e6 * 0.01,
    local uB = 2 * F / EA * 1 / 3 * l,
    local uC = F / EA * 2 / 5 * l,

    primary: {
      B: test.Ux(if hinge == 'left' then -uB else uC),
    },

    interpolation: {
      AB: test.Constant('Nx', 2 * F, range=[0, l / 3]) +
          test.Constant('Nx', 0, range=[l / 3, l]) +
          test.Linear('Ux', 0, uB, range=[0, l / 3]) +
          test.Constant('Ux', uB, range=[l / 3, l]),
      BC: test.Constant('Nx', 0, range=[0, 3 / 5 * l]) +
          test.Constant('Nx', F, range=[3 / 5 * l, l]) +
          test.Constant('Ux', -uC, range=[0, 3 / 5 * l]) +
          test.Linear('Ux', -uC, 0, range=[3 / 5 * l, l]),
    },
  },
};

[
  horizontal_axial_hinge(hinge='left'),
  horizontal_axial_hinge(hinge='right'),
]
