local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local common(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [2 * l, 0, 0],
    D: [3 * l, 0, 0],
  },

  material: bvp.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=10e-6),

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Phiy'], C: ['Phiy'] }),
    CD: bvp.Frame2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Uz() + bvp.Phiy(),
  },
};

local phi_phi_fz(F, l) = common(l) {
  name: 'phi_phi_fz_load_attribution_%g_%g' % [F, l],

  neumann: {
    AB: test.Fz(F / 3, l),
    BC: test.Fz(F / 3, 0) + test.Fz(F / 3, l),
    CD: test.Fz(F / 3, 0),
    B: test.Fz(-F / 3),
    C: test.Fz(-F / 3),
  },

  expected: {
    reaction: {
      A: test.Fz(F) + test.Fx(0) + test.My(-F * l),
      D: test.Fz(F) + test.My(F * l),
    },
    interpolation: {
      AB: test.Constant('Vz', F) + test.Linear('My', -F * l, 0),
      BC: test.Constant('Vz', 0) + test.Constant('My', 0),
      CD: test.Constant('Vz', -F) + test.Linear('My', 0, -F * l),
    },
  },
};

local phi_phi_my(M, l) = common(l) {
  name: 'phi_phi_my_load_attribution_%g_%g' % [M, l],

  neumann: {
    AB: bvp.My(M / 3, l),
    BC: bvp.My(M / 3, 0) + bvp.My(M / 3, l),
    CD: bvp.My(M / 3, 0),
    B: bvp.My(-M / 3),
    C: bvp.My(-M / 3),
  },

  expected: {
    reaction: {
      A: test.Fz(2 / 3 * M / l) + test.Fx(0) + test.My(0),
      D: test.Fz(-2 / 3 * M / l) + test.My(0),
    },
    interpolation: {
      AB: test.Constant('Vz', 2 / 3 * M / l) + test.Linear('My', 0, 2 / 3 * M),
      BC: test.Constant('Vz', 2 / 3 * M / l) + test.Linear('My', -M / 3, M / 3),
      CD: test.Constant('Vz', 2 / 3 * M / l) + test.Linear('My', -2 / 3 * M, 0),
    },
  },
};

[
  phi_phi_fz(F=5e3, l=5),
  phi_phi_my(M=1.2e3, l=4),
]
