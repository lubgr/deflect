local lib = import 'common.libsonnet';

local default = lib.Defaults();

local common(l, E, Iyy) = {
  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  material: lib.LinElast('default', E=E, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=1, Iyy=Iyy),

  elements: {
    AB: default.Frame2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz() + lib.Phiy(),
  },
};

local with_nodal_fz(F, l, E, Iyy) = common(l, E, Iyy) {
  name: 'nodal_fz_%g' % F,

  neumann: {
    B: lib.Fz(F),
  },

  expected: {
    reaction: {
      A: lib.Fx(0) + lib.Fz(-F),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(F * std.pow(l, 3) / (3 * E * Iyy)) +
         lib.Phiy(-F * l * l / (2 * E * Iyy)),
    },
  },
};

local with_nodal_my(M, l, E, Iyy) = common(l, E, Iyy) {
  name: 'nodal_my_%g' % M,

  neumann: {
    B: lib.My(M),
  },

  expected: {
    reaction: {
      A: lib.Fx(0) + lib.My(-M),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(-M * std.pow(l, 2) / (2 * E * Iyy)) +
         lib.Phiy(M * l / (E * Iyy)),
    },
  },
};

local with_constant(q, l, E, Iyy) = common(l, E, Iyy) {
  name: 'constant_qz_%g' % q,

  neumann: {
    AB: lib.qz(q),
  },

  expected: {
    reaction: {
      A: lib.Fx(0) + lib.Fz(q * l) + lib.My(-q * l * l / 2),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(-q * std.pow(l, 4) / (8 * E * Iyy)) +
         lib.Phiy(q * std.pow(l, 3) / (6 * E * Iyy)),
    },
  },
};

[
  with_nodal_fz(F=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_fz(F=1e3, l=5.75, E=10000e6, Iyy=10e-6),

  with_nodal_my(M=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_my(M=1e3, l=5.75, E=10000e6, Iyy=10e-6),

  with_constant(q=-1e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_constant(q=1e3, l=5.75, E=10000e6, Iyy=10e-6),
]
