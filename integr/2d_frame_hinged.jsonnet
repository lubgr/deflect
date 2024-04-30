local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local frame(l, h, Ih, Iv) = {
  nodes:
    {
      A: [0, 0, 0],
      B: [l, 0, 0],
      C: [0, 0, h],
      D: [l, 0, h],
    },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
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
    A: bvp.Ux() + bvp.Uz(),
    B: bvp.Ux() + bvp.Uz(),
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
    TOP: bvp.qz(q),
  },

  expected: relaxed_tol {
    local H = q * l * l / (4 * h * (2 * $.k + 3)),

    reaction: {
      A: test.Fz(q * l / 2) + test.Fx(H),
      B: test.Fz(q * l / 2) + test.Fx(-H),
    },
    interpolation: {
      LEFT: test.Linear('My', 0, -H * h),
      TOP: test.Quadratic('My', eval=[[0, -H * h], [l / 2, -H * h + q * l * l / 8], [l, -H * h]]),
      RIGHT: test.Linear('My', -H * h, 0),
    },
  },
};

local const_qz_right(q, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'const_qz_right_%g' % q,

  neumann: {
    RIGHT: bvp.qz(q),
  },

  expected: relaxed_tol {
    local H1 = q * h / 8 * (5 * $.k + 6) / (2 * $.k + 3),
    local H2 = H1 - q * h,

    reaction: {
      local V = q * h * h / (2 * l),
      A: test.Fz(V) + test.Fx(H1),
      B: test.Fz(-V) + test.Fx(-H2),
    },
    interpolation: {
      LEFT: test.Linear('My', 0, -H1 * h),
      TOP: test.Linear('My', -H1 * h, -H2 * h - q * h * h / 2),
      RIGHT: test.Quadratic('My', eval=[[0, -H2 * h - q * h * h / 2], [h, 0]]),
    },
  },
};

local fz_element_top(F, a, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'fz_top_%g_%g' % [F, a],

  neumann: {
    TOP: bvp.Fz(F, a),
  },

  expected: relaxed_tol {
    local b = l - a,
    local H = 3 / 2 * F * a * b / (h * l * (2 * $.k + 3)),

    reaction: {
      A: test.Fz(F * b / l) + test.Fx(H),
      B: test.Fz(F * a / l) + test.Fx(-H),
    },
    interpolation: {
      LEFT: test.Linear('My', 0, -H * h),
      TOP: test.Linear('My', -H * h, -H * h + a * b * F / l, range=[0, a]) +
           test.Linear('My', -H * h + a * b * F / l, -H * h, range=[a, l]),
      RIGHT: test.Linear('My', -H * h, 0),
    },
  },
};

local fx_node_top_right(F, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  name: 'fz_node_right_top_%g' % F,

  neumann: {
    D: bvp.Fx(-F),
  },

  expected: relaxed_tol {
    reaction: {
      A: test.Fx(F / 2) + test.Fz(F * h / l),
      B: test.Fx(F / 2) + test.Fz(-F * h / l),
    },
    interpolation: {
      local M = -F / 2 * h,
      LEFT: test.Linear('My', 0, M),
      TOP: test.Linear('My', M, -M),
      RIGHT: test.Linear('My', -M, 0),
    },
  },
};

local fx_right(F, b, l, h, Ih, Iv) = frame(l, h, Ih, Iv) + {
  local a = h - b,
  name: 'fz_right_%g_%g' % [F, a],

  neumann: {
    // Need to convert between a and b since the element goes from top to bottom:
    RIGHT: test.Fz(F, b),
  },

  expected: relaxed_tol {
    local H1 = 3 * F * a * $.k / (2 * h * (2 * $.k + 3)) * (1 - a * a / (3 * h * h) + 1 / $.k),
    local H2 = H1 - F,

    reaction: {
      A: test.Fz(F * a / l) + test.Fx(H1),
      B: test.Fz(-F * a / l) + test.Fx(-H2),
    },
    interpolation: {
      LEFT: test.Linear('My', 0, -H1 * h),
      TOP: test.Linear('My', -H1 * h, -H2 * h - F * (h - a)),
      RIGHT: test.Linear('My', -H2 * h - F * (h - a), -H2 * a, range=[0, b]) +
             test.Linear('My', -H2 * a, 0, range=[b, h]),
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
