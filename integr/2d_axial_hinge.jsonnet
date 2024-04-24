local lib = import 'common.libsonnet';

local default = lib.Defaults();

local with_axial_hinge(q, F, l) = {
  name: 'with_axial_hinge',
  description: '3 frames with one normal force hinge',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, -l],
    C: [2 * l, 0, 0],
    D: [l, 0, 0],
  },

  material: lib.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=0.01, Iyy=10e-6),

  elements: {
    AD: default.Frame2d(),
    DC: default.Frame2d(hinges={ D: ['Ux'] }),
    DB: default.Frame2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Uz(),
    C: lib.Ux(),
  },

  neumann: {
    B: lib.Fx(-F),
    C: lib.Fz(-F),
    AD: lib.qx(-q),
    DC: lib.qx(-q),
  },

  expected: {
    reaction: {
      A: lib.Fx(F + q * l) + lib.Fz(-2 * F),
      B: lib.Fz(3 * F),
      C: lib.Fx(q * l),
    },
    interpolation: {
      AD: lib.Linear('Nx', -F - q * l, -F) +
          lib.Constant('Vz', -2 * F) +
          lib.Linear('My', 0, -2 * F * l),
      DC: lib.Linear('Nx', 0, q * l) +
          lib.Constant('Vz', F) +
          lib.Linear('My', -F * l, 0),
      DB: lib.Constant('Nx', -3 * F) +
          lib.Constant('Vz', F) +
          lib.Linear('My', -F * l, 0),
    },
  },
};

[
  with_axial_hinge(2e3, 5e3, 2.5),
  with_axial_hinge(-1e3, -2e3, 0.5),
]
