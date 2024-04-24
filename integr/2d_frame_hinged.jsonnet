local lib = import 'common.libsonnet';

local frame(l, h, Ih, Iv) = {
  nodes:
    {
      A: [0, 0, 0],
      B: [l, 0, 0],
      C: [0, 0, h],
      D: [l, 0, h],
    },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: {
    beam: {
      kind: 'constants',
      parameter: {
        A: 0.01,
        Iyy: Ih,
      },
    },
    pillar: {
      kind: 'constants',
      parameter: {
        A: 0.01,
        Iyy: Iv,
      },
    },
  },

  local frame(n0, n1, cs) = {
    kind: 'frame2d',
    material: 'default',
    cs: cs,
    nodes: [n0, n1],
  },

  elements: {
    LEFT: frame('A', 'C', 'pillar'),
    RIGHT: frame('D', 'B', 'pillar'),
    TOP: frame('C', 'D', 'beam'),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Ux() + lib.Uz(),
  },

  // Constant that helps expressing expected reactions:
  k:: Ih / Iv * h / l,

};

local relaxed_tol = {
  tolerance: {
    reaction: 5e-3,
    polynomial: 5e-3,
  },
};

local const_qz_top(q, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'const_qz_top_%g' % q,

  neumann: {
    TOP: lib.qz(q),
  },

  expected: relaxed_tol {
    local H = q * l * l / (4 * h * (2 * $.k + 3)),

    reaction: {
      A: lib.Fz(q * l / 2) + lib.Fx(H),
      B: lib.Fz(q * l / 2) + lib.Fx(-H),
    },
    interpolation: {
      LEFT: lib.Linear('My', 0, -H * h),
      TOP: lib.Quadratic('My', eval=[[0, -H * h], [l / 2, -H * h + q * l * l / 8], [l, -H * h]]),
      RIGHT: lib.Linear('My', -H * h, 0),
    },
  },
};

local const_qz_right(q, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'const_qz_right_%g' % q,

  neumann: {
    RIGHT: lib.qz(q),
  },

  expected: relaxed_tol {
    local H1 = q * h / 8 * (5 * $.k + 6) / (2 * $.k + 3),
    local H2 = H1 - q * h,

    reaction: {
      local V = q * h * h / (2 * l),
      A: lib.Fz(V) + lib.Fx(H1),
      B: lib.Fz(-V) + lib.Fx(-H2),
    },
    interpolation: {
      LEFT: lib.Linear('My', 0, -H1 * h),
      TOP: lib.Linear('My', -H1 * h, -H2 * h - q * h * h / 2),
      RIGHT: lib.Quadratic('My', eval=[[0, -H2 * h - q * h * h / 2], [h, 0]]),
    },
  },
};

local fz_element_top(F, a, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'fz_top_%g_%g' % [F, a],

  neumann: {
    TOP: lib.Fz(F, a),
  },

  expected: relaxed_tol {
    local b = l - a,
    local H = 3 / 2 * F * a * b / (h * l * (2 * $.k + 3)),

    reaction: {
      A: lib.Fz(F * b / l) + lib.Fx(H),
      B: lib.Fz(F * a / l) + lib.Fx(-H),
    },
    interpolation: {
      LEFT: lib.Linear('My', 0, -H * h),
      TOP: lib.Linear('My', -H * h, -H * h + a * b * F / l, range=[0, a]) +
           lib.Linear('My', -H * h + a * b * F / l, -H * h, range=[a, l]),
      RIGHT: lib.Linear('My', -H * h, 0),
    },
  },
};

local fx_node_top_right(F, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'fz_node_right_top_%g' % F,

  neumann: {
    D: lib.Fx(-F),
  },

  expected: relaxed_tol {
    reaction: {
      A: lib.Fx(F / 2) + lib.Fz(F * h / l),
      B: lib.Fx(F / 2) + lib.Fz(-F * h / l),
    },
    interpolation: {
      local M = -F / 2 * h,
      LEFT: lib.Linear('My', 0, M),
      TOP: lib.Linear('My', M, -M),
      RIGHT: lib.Linear('My', -M, 0),
    },
  },
};

local fx_right(F, b, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  local a = h - b,
  name: 'fz_right_%g_%g' % [F, a],

  neumann: {
    // Need to convert between a and b since the element goes from top to bottom:
    RIGHT: lib.Fz(F, b),
  },

  expected: relaxed_tol {
    local H1 = 3 * F * a * $.k / (2 * h * (2 * $.k + 3)) * (1 - a * a / (3 * h * h) + 1 / $.k),
    local H2 = H1 - F,

    reaction: {
      A: lib.Fz(F * a / l) + lib.Fx(H1),
      B: lib.Fz(-F * a / l) + lib.Fx(-H2),
    },
    interpolation: {
      LEFT: lib.Linear('My', 0, -H1 * h),
      TOP: lib.Linear('My', -H1 * h, -H2 * h - F * (h - a)),
      RIGHT: lib.Linear('My', -H2 * h - F * (h - a), -H2 * a, range=[0, b]) +
             lib.Linear('My', -H2 * a, 0, range=[b, h]),
    },
  },
};


[
  const_qz_top(q=10e3, l=3.5, h=1.8, Ih=10e-6, Iv=20e-6),
  const_qz_right(q=10e3, l=3.5, h=1.8, Ih=10e-6, Iv=20e-6),
  fz_element_top(F=50e3, a=1.25, l=3, h=2, Ih=10e-6, Iv=20e-6),
  fz_element_top(F=25e3, a=2.9, l=3, h=3.75, Ih=20e-6, Iv=10e-6),
  fx_node_top_right(F=25e3, l=2, h=5, Ih=10e-6, Iv=12.5e-6),
  fx_right(F=50e3, b=1.25, l=4, h=2.5, Ih=10e-6, Iv=20e-6),
  fx_right(F=20e3, b=0.5, l=2, h=2.5, Ih=20e-6, Iv=15e-6),
]
