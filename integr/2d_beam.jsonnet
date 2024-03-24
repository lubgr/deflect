local lib = import 'common.libsonnet';

local default = lib.Defaults();

local common(l, E, Iyy) = {
  material: lib.LinElast('default', E=E, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=1, Iyy=Iyy),

  elements: {
    AB: default.Frame2d(),
  },
};

local horizontal(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Uz(),
  },
};

local vertical(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [0, 0, l],
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Ux(),
  },
};

local with_const_qz(q, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) {
  name: 'element_const_qz_%g' % q,

  neumann: {
    AB: lib.qz(q),
  },

  expected: {
    primary: {
      A: lib.Phiy(q * std.pow(l, 3) / (24 * E * Iyy)),
      B: lib.Phiy(-q * std.pow(l, 3) / (24 * E * Iyy)),
    },
    reaction: {
      A: lib.Fz(q * l / 2),
      B: lib.Fz(q * l / 2),
    },
  },
};

local with_const_qz_vertical(q, l, E, Iyy) = with_const_qz(q, l, E, Iyy) + vertical(l) + {
  name: 'element_const_qz_%g_vertical' % q,

  expected: super.expected + {
    primary: super.primary,
    reaction: {
      A: lib.Fx(-q * l / 2),
      B: lib.Fx(-q * l / 2),
    },
  },
};

local with_linear_qz(q0, q1, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'element_linear_qz_%.1f_%.1f' % [q0, q1],

  neumann: {
    AB: lib.qz([q0, q1]),
  },

  expected: {
    primary: {
      A: lib.Phiy((8 * q0 + 7 * q1) * std.pow(l, 3) / (360 * E * Iyy)),
      B: lib.Phiy(-(7 * q0 + 8 * q1) * std.pow(l, 3) / (360 * E * Iyy)),
    },
    reaction: {
      A: lib.Fz((2 * q0 + q1) * l / 6),
      B: lib.Fz((q0 + 2 * q1) * l / 6),
    },
  },
};

local with_element_fz(F, x, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'element_fz_%.1f_%.2f' % [F, x],

  neumann: {
    AB: lib.Fz(F, x),
  },

  expected: {
    primary: {
      A: lib.Phiy(F * x * (l - x) * (2 * l - x) / (6 * E * Iyy * l)),
      B: lib.Phiy(-F * x * (l - x) * (l + x) / (6 * E * Iyy * l)),
    },
    reaction: {
      A: lib.Fz(F * (1 - x / l)),
      B: lib.Fz(F * x / l),
    },
  },
};

local with_element_fz_vertical(F, x, l, E, Iyy) = with_element_fz(F, x, l, E, Iyy) + vertical(l) + {
  name: 'element_fz_%g_%.2f_vertical' % [F, x],

  expected: super.expected + {
    primary: super.primary,
    reaction: {
      A: lib.Fx(-F * (1 - x / l)),
      B: lib.Fx(-F * x / l),
    },
  },
};

local with_all(l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'all_loads_01',

  local fconcentrated = with_element_fz(11.5e3, 0.3 * l, l, E, Iyy),
  local qconst = with_const_qz(1.25e3, l, E, Iyy),
  local qlinear = with_linear_qz(-1.25e3, 1.25e3, l, E, Iyy),

  neumann: {
    AB: fconcentrated.neumann.AB + qconst.neumann.AB + qlinear.neumann.AB,
  },

  expected: {
    primary: {
      local phiy(from, node) = from.expected.primary[node][0].Phiy,
      A: lib.Ux(0) + lib.Uz(0) +
         lib.Phiy(phiy(fconcentrated, 'A') + phiy(qconst, 'A') + phiy(qlinear, 'A')),
      B: lib.Ux(0) + lib.Uz(0) +
         lib.Phiy(phiy(fconcentrated, 'B') + phiy(qconst, 'B') + phiy(qlinear, 'B')),
    },
    reaction: {
      local fz(from, node) = from.expected.reaction[node][0].Fz,
      A: lib.Fx(0) +
         lib.Fz(fz(fconcentrated, 'A') + fz(qconst, 'A') + fz(qlinear, 'A')),
      B: lib.Fz(fz(fconcentrated, 'B') + fz(qconst, 'B') + fz(qlinear, 'B')),
    },
  },
};

[
  with_const_qz(q=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_const_qz(q=2.5e3, l=6, E=12345e6, Iyy=8e-6),
  with_const_qz_vertical(q=2.5e3, l=6, E=12345e6, Iyy=8e-6),

  with_linear_qz(q0=1e3, q1=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=0e3, q1=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=1e3, q1=0e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=-2.5e3, q1=10e3, l=2, E=12345e6, Iyy=8e-6),
  with_linear_qz(q0=2.5e3, q1=-8e3, l=2, E=12345e6, Iyy=8e-6),
]
+
[
  with_element_fz(F=10e3, x=x, l=4, E=10000e6, Iyy=10e-6)
  for x in [0, 0.25, 2, 3.9, 4]
]
+
[
  with_element_fz_vertical(F=10e3, x=x, l=4, E=10000e6, Iyy=10e-6)
  for x in [0, 0.25, 2, 3.9, 4]
]
+
[
  with_all(l=4, E=10000e6, Iyy=10e-6),
]
