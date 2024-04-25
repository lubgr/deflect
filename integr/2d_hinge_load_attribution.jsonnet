local lib = import 'common.libsonnet';

local default = lib.Defaults();

local common(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [2 * l, 0, 0],
    D: [3 * l, 0, 0],
  },

  material: lib.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=0.01, Iyy=10e-6),

  elements: {
    AB: default.Frame2d(),
    BC: default.Frame2d(hinges={ B: ['Phiy'], C: ['Phiy'] }),
    CD: default.Frame2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz() + lib.Phiy(),
    D: lib.Uz() + lib.Phiy(),
  },
};

local phi_phi_fz(F, l) = common(l) {
  name: 'phi_phi_fz_load_attribution_%g_%g' % [F, l],

  neumann: {
    AB: lib.Fz(F / 3, l),
    BC: lib.Fz(F / 3, 0) + lib.Fz(F / 3, l),
    CD: lib.Fz(F / 3, 0),
    B: lib.Fz(-F / 3),
    C: lib.Fz(-F / 3),
  },

  expected: {
    reaction: {
      A: lib.Fz(F) + lib.Fx(0) + lib.My(-F * l),
      D: lib.Fz(F) + lib.My(F * l),
    },
    interpolation: {
      AB: lib.Constant('Vz', F) + lib.Linear('My', -F * l, 0),
      BC: lib.Constant('Vz', 0) + lib.Constant('My', 0),
      CD: lib.Constant('Vz', -F) + lib.Linear('My', 0, -F * l),
    },
  },
};

local phi_phi_my(M, l) = common(l) {
  name: 'phi_phi_my_load_attribution_%g_%g' % [M, l],

  neumann: {
    AB: lib.My(M / 3, l),
    BC: lib.My(M / 3, 0) + lib.My(M / 3, l),
    CD: lib.My(M / 3, 0),
    B: lib.My(-M / 3),
    C: lib.My(-M / 3),
  },

  expected: {
    reaction: {
      A: lib.Fz(2 / 3 * M / l) + lib.Fx(0) + lib.My(0),
      D: lib.Fz(-2 / 3 * M / l) + lib.My(0),
    },
    interpolation: {
      AB: lib.Constant('Vz', 2 / 3 * M / l) + lib.Linear('My', 0, 2 / 3 * M),
      BC: lib.Constant('Vz', 2 / 3 * M / l) + lib.Linear('My', -M / 3, M / 3),
      CD: lib.Constant('Vz', 2 / 3 * M / l) + lib.Linear('My', -2 / 3 * M, 0),
    },
  },
};

[
  phi_phi_fz(F=5e3, l=5),
  phi_phi_my(M=1.2e3, l=4),
]
