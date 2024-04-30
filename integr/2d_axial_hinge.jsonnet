local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

local with_axial_hinge(q, F, l) = {
  name: 'with_axial_hinge',
  description: '3 frames with one normal force hinge',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, -l],
    C: [2 * l, 0, 0],
    D: [l, 0, 0],
  },

  material: bvp.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=10e-6),

  elements: {
    AD: default.Frame2d(),
    DC: default.Frame2d(hinges={ D: ['Ux'] }),
    DB: default.Frame2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    B: bvp.Uz(),
    C: bvp.Ux(),
  },

  neumann: {
    B: bvp.Fx(-F),
    C: bvp.Fz(-F),
    AD: bvp.qx(-q),
    DC: bvp.qx(-q),
  },

  expected: {
    reaction: {
      A: test.Fx(F + q * l) + test.Fz(-2 * F),
      B: test.Fz(3 * F),
      C: test.Fx(q * l),
    },
    interpolation: {
      AD: test.Linear('Nx', -F - q * l, -F) +
          test.Constant('Vz', -2 * F) +
          test.Linear('My', 0, -2 * F * l),
      DC: test.Linear('Nx', 0, q * l) +
          test.Constant('Vz', F) +
          test.Linear('My', -F * l, 0),
      DB: test.Constant('Nx', -3 * F) +
          test.Constant('Vz', F) +
          test.Linear('My', -F * l, 0),
    },
  },
};

[
  with_axial_hinge(2e3, 5e3, 2.5),
  with_axial_hinge(-1e3, -2e3, 0.5),
]
