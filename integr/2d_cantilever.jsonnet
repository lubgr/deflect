local lib = import 'common.libsonnet';

local default = lib.Defaults();

local common(l, E, Iyy, fix) = {
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
    assert fix == 'left' || fix == 'right',
    [if fix == 'left' then 'A' else 'B']: lib.Ux() + lib.Uz() + lib.Phiy(),
  },
};

local with_nodal_fz(F, l, E, Iyy) = common(l, E, Iyy, fix='left') {
  name: 'nodal_fz_%g' % F,

  neumann: {
    B: lib.Fz(F),
  },

  expected: {
    local phiy(x) = F / (E * Iyy) * (l * x - x * x / 2),
    local uz(x) = -F / (E * Iyy) * (l / 2 * std.pow(x, 2) - std.pow(x, 3) / 6),

    reaction: {
      A: lib.Fx(0) + lib.Fz(-F),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(-uz(l)) +
         lib.Phiy(-phiy(l)),
    },
    interpolation: {
      AB: lib.Constant('Vz', -F) +
          lib.Linear('My', F * l, 0) +
          lib.Quadratic('Phiy', eval=lib.Samples(phiy, 0, l, 5)) +
          lib.Cubic('Uz', eval=lib.Samples(uz, 0, l, 5)),
    },
  },
};

local with_nodal_fz_inverse(F, l, E, Iyy) = common(l, E, Iyy, fix='right') {
  name: 'nodal_fz_%g_inverse' % F,

  neumann: {
    A: lib.Fz(F),
  },

  expected: {
    local phiy(x) = F / (E * Iyy) * (x * x / 2 - l * l / 2),
    local uz(x) = -F / (E * Iyy) * (std.pow(x, 3) / 6 - l * l * x / 2 + std.pow(l, 3) / 3),

    reaction: {
      B: lib.Fx(0) + lib.Fz(-F),
    },
    primary: {
      A: lib.Uz(-uz(0)) + lib.Phiy(-phiy(0)),
      B: lib.Uz(-uz(l)) + lib.Phiy(-phiy(l)),
    },
    interpolation: {
      AB: lib.Constant('Vz', F) +
          lib.Linear('My', 0, F * l) +
          lib.Quadratic('Phiy', eval=lib.Samples(phiy, 0, l, 5)) +
          lib.Cubic('Uz', eval=lib.Samples(uz, 0, l, 5)),
    },
  },
};

local with_nodal_my(M, l, E, Iyy) = common(l, E, Iyy, fix='left') {
  name: 'nodal_my_%g' % M,

  neumann: {
    B: lib.My(M),
  },

  expected: {
    local uz(x) = M / (E * Iyy) * (x * x / 2),

    reaction: {
      A: lib.Fz(0) + lib.My(-M),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(-uz(l)) +
         lib.Phiy(M * l / (E * Iyy)),
    },
    interpolation: {
      AB: lib.Constant('Vz', 0) +
          lib.Constant('My', -M) +
          lib.Quadratic('Uz', eval=lib.Samples(uz, 0, l, 5)),
    },
  },
};

local with_nodal_my_inverse(M, l, E, Iyy) = common(l, E, Iyy, fix='right') {
  name: 'nodal_my_inverse_%g' % M,

  neumann: {
    A: lib.My(M),
  },

  expected: {
    local uz(x) = -M / (E * Iyy) * (l * l / 2 - l * x + x * x / 2),

    reaction: {
      B: lib.Fz(0) + lib.My(-M),
    },
    primary: {
      A: lib.Uz(-uz(0)),
    },
    interpolation: {
      AB: lib.Constant('Vz', 0) +
          lib.Constant('My', M) +
          lib.Quadratic('Uz', eval=lib.Samples(uz, 0, l, 5)),
    },
  },
};

local with_constant(q, l, E, Iyy) = common(l, E, Iyy, fix='left') {
  name: 'constant_qz_%g' % q,

  neumann: {
    AB: lib.qz(q),
  },

  expected: {
    local uz(x) = q / (E * Iyy) * (std.pow(x, 4) / 24 - std.pow(x, 3) * l / 6 + std.pow(l * x, 2) / 4),

    reaction: {
      A: lib.Fx(0) + lib.Fz(q * l) + lib.My(-q * l * l / 2),
    },
    primary: {
      B: lib.Ux(0) +
         lib.Uz(-uz(l)) +
         lib.Phiy(q * std.pow(l, 3) / (6 * E * Iyy)),
    },
    interpolation: {
      AB: lib.Linear('Vz', q * l, 0) +
          lib.Quadratic('My', eval=lib.Samples(function(x) (-q / 2 * std.pow(l - x, 2)), 0, l, 5)) +
          lib.Quartic('Uz', eval=lib.Samples(uz, 0, l, 7)),
    },
  },
};

local with_constant_inverse(q, l, E, Iyy) = common(l, E, Iyy, fix='right') {
  name: 'constant_qz_%g_inverse' % q,

  neumann: {
    AB: lib.qz(q),
  },

  expected: {
    local uz(x) = q * std.pow(l, 4) / (24 * E * Iyy) * (3 - 4 * x / l + std.pow(x / l, 4)),

    reaction: {
      B: lib.Fx(0) + lib.Fz(q * l) + lib.My(q * l * l / 2),
    },
    primary: {
      A: lib.Uz(-uz(0)),
    },
    interpolation: {
      AB: lib.Linear('Vz', 0, -q * l) +
          lib.Quadratic('My', eval=[[0, 0], [l, -q * l * l / 2]]) +
          lib.Quartic('Uz', eval=lib.Samples(uz, 0, l, 7)),
    },
  },
};

[
  with_nodal_fz(F=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_fz(F=1e3, l=5.75, E=10000e6, Iyy=10e-6),
  with_nodal_fz_inverse(F=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_fz_inverse(F=1e3, l=5.75, E=10000e6, Iyy=10e-6),

  with_nodal_my(M=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_my(M=1e3, l=5.75, E=10000e6, Iyy=10e-6),
  with_nodal_my_inverse(M=-2.5e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_nodal_my_inverse(M=1e3, l=5.75, E=10000e6, Iyy=10e-6),

  with_constant(q=-1e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_constant(q=1e3, l=5.75, E=10000e6, Iyy=10e-6),
  with_constant_inverse(q=-1e3, l=1, E=12345e6, Iyy=98.7655e-6),
  with_constant_inverse(q=1e3, l=5.75, E=10000e6, Iyy=10e-6),
]
