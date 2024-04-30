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
};

local qz_common(q, l) = common(l) {
  local elements = self.elements,

  neumann: {
    [id]: bvp.qz(q)
    for id in std.objectFields(elements)
  },
};

local qz_hinge_phi_rightward(q, l) = qz_common(q, l) + {
  name: 'qz_hinge_phi_rightward',

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Uz(),
  },

  elements: {
    AB: bvp.Frame2d(hinges={ B: ['Phiy'] }),
    BC: bvp.Frame2d(),
    CD: bvp.Frame2d(),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(2 * q * l) + test.My(-3 / 2 * q * l * l),
      D: test.Fz(q * l),
    },
    interpolation: {
      AB: test.Quadratic('My', eval=[[0, -3 / 2 * q * l * l], [l, 0]]),
      BC: test.Quadratic('My', eval=[[0, 0], [l, q * l * l / 2]]),
      CD: test.Quadratic('My', eval=[[0, q * l * l / 2], [l, 0]]),
    },
  },
};

local qz_hinge_phi_leftward(q, l) = qz_hinge_phi_rightward(q, l) + qz_common(-q, l) + {
  name: 'qz_hinge_phi_leftward',

  elements: {
    BA: bvp.Frame2d(hinges={ B: ['Phiy'] }),
    CB: bvp.Frame2d(),
    DC: bvp.Frame2d(),
  },

  expected: super.expected {
    interpolation: {
      BA: test.Quadratic('My', eval=[[0, 0], [l, 3 / 2 * q * l * l]]),
      CB: test.Quadratic('My', eval=[[0, -q * l * l / 2], [l, 0]]),
      DC: test.Quadratic('My', eval=[[0, 0], [l, -q * l * l / 2]]),
    },
  },
};

local qz_hinge_uz_rightward(q, l) = qz_common(q, l) + {
  name: 'qz_hinge_uz_rightward',

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Uz(),
  },

  elements: {
    AB: bvp.Frame2d(hinges={ B: ['Uz'] }),
    BC: bvp.Frame2d(),
    CD: bvp.Frame2d(),
  },

  expected: {
    local l2 = l * l,

    reaction: {
      A: test.Fx(0) + test.Fz(q * l) + test.My(3 / 2 * q * l2),
      D: test.Fz(2 * q * l),
    },
    interpolation: {
      AB: test.Quadratic('My', eval=[[0, 3 / 2 * q * l2], [l, 2 * q * l2]]),
      BC: test.Quadratic('My', eval=[[0, 2 * q * l2], [l, 3 / 2 * q * l2]]),
      CD: test.Quadratic('My', eval=[[0, 3 / 2 * q * l2], [l, 0]]),
    },
  },
};

local qz_hinge_uz_leftward(q, l) = qz_hinge_uz_rightward(q, l) + qz_common(-q, l) + {
  name: 'qz_hinge_uz_leftward',

  elements: {
    BA: bvp.Frame2d(hinges={ B: ['Uz'] }),
    CB: bvp.Frame2d(),
    DC: bvp.Frame2d(),
  },

  expected: super.expected {
    local l2 = l * l,

    interpolation: {
      BA: test.Quadratic('My', eval=[[0, -2 * q * l2], [l, -3 / 2 * q * l2]]),
      CB: test.Quadratic('My', eval=[[0, -3 / 2 * q * l2], [l, -2 * q * l2]]),
      DC: test.Quadratic('My', eval=[[0, 0], [l, -3 / 2 * q * l2]]),
    },
  },
};

local qz_hinge_phi_phi(q, l) = qz_common(q, l) + {
  name: 'qz_hinge_phi_phi',

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Uz() + bvp.Phiy(),
  },

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Phiy'], C: ['Phiy'] }),
    CD: bvp.Frame2d(),
  },

  expected: {
    local l2 = l * l,

    reaction: {
      A: test.Fx(0) + test.Fz(3 / 2 * q * l) + test.My(-q * l2),
      D: test.Fz(3 / 2 * q * l) + test.My(q * l2),
    },
    interpolation: {
      AB: test.Quadratic('My', eval=[[0, -q * l2], [l, 0]]),
      BC: test.Quadratic('My', eval=[[0, 0], [l / 2, q * l2 / 8], [l, 0]]),
      CD: test.Quadratic('My', eval=[[0, 0], [l, -q * l2]]),
    },
  },
};

local qz_hinge_phi_phi_variation_1(q, l) = qz_hinge_phi_phi(q, l) + {
  name: 'qz_hinge_phi_phi_variation_1',

  elements: {
    AB: bvp.Frame2d(hinges={ B: ['Phiy'] }),
    BC: bvp.Frame2d(hinges={ C: ['Phiy'] }),
    CD: bvp.Frame2d(),
  },
};

local qz_hinge_phi_phi_variation_2(q, l) = qz_hinge_phi_phi(q, l) + {
  name: 'qz_hinge_phi_phi_variation_2',

  elements: {
    AB: bvp.Frame2d(hinges={ B: ['Phiy'] }),
    BC: bvp.Frame2d(),
    CD: bvp.Frame2d(hinges={ C: ['Phiy'] }),
  },
};

local qz_hinge_phi_phi_variation_3(q, l) = qz_hinge_phi_phi(q, l) + {
  name: 'qz_hinge_phi_phi_variation_3',

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Phiy'] }),
    CD: bvp.Frame2d(hinges={ C: ['Phiy'] }),
  },
};

local qz_hinge_uz_phi_rightward(q, l) = qz_common(q, l) {
  name: 'qz_hinge_uz_phi_rightward',

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Uz() + bvp.Phiy(),
  },

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Uz'], C: ['Phiy'] }),
    CD: bvp.Frame2d(),
  },

  expected: {
    local l2 = l * l,

    reaction: {
      A: test.Fx(0) + test.Fz(q * l) + test.My(0),
      D: test.Fz(2 * q * l) + test.My(3 / 2 * q * l2),
    },
    interpolation: {
      AB: test.Quadratic('My', eval=[[0, 0], [l, q * l2 / 2]]),
      BC: test.Quadratic('My', eval=[[0, q * l2 / 2], [l, 0]]),
      CD: test.Quadratic('My', eval=[[0, 0], [l, -3 / 2 * q * l2]]),
    },
  },
};

local qz_hinge_uz_phi_leftward(q, l) = qz_hinge_uz_phi_rightward(q, l) + qz_common(-q, l) {
  name: 'qz_hinge_uz_phi_leftward',

  elements: {
    BA: bvp.Frame2d(),
    CB: bvp.Frame2d(hinges={ B: ['Uz'], C: ['Phiy'] }),
    DC: bvp.Frame2d(),
  },

  expected: super.expected {
    local l2 = l * l,

    interpolation: {
      BA: test.Quadratic('My', eval=[[0, -q * l2 / 2], [l, 0]]),
      CB: test.Quadratic('My', eval=[[0, 0], [l, -q * l2 / 2]]),
      DC: test.Quadratic('My', eval=[[0, 3 / 2 * q * l2], [l, 0]]),
    },
  },
};

local fz_common(F, l) = common(l) {
  neumann: {
    B: bvp.Fz(-F),
    C: bvp.Fz(-F),
  },
};

local fz_hinge_phi_rightward(F, l) = fz_common(F, l) {
  name: 'fz_hinge_phi_rightward',
  dirichlet: qz_hinge_phi_rightward(0, l).dirichlet,
  elements: qz_hinge_phi_rightward(0, l).elements,

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(3 / 2 * F) + test.My(-3 / 2 * F * l),
      D: test.Fz(F / 2),
    },
    interpolation: {
      AB: test.Linear('My', -3 / 2 * F * l, 0),
      BC: test.Linear('My', 0, F * l / 2),
      CD: test.Linear('My', F * l / 2, 0),
    },
  },
};

local fz_hinge_phi_leftward(F, l) = fz_hinge_phi_rightward(F, l) {
  name: 'fz_hinge_phi_leftward',
  elements: qz_hinge_phi_leftward(0, l).elements,

  expected: super.expected + {
    interpolation: {
      BA: test.Linear('My', 0, 3 / 2 * F * l),
      CB: test.Linear('My', -F * l / 2, 0),
      DC: test.Linear('My', 0, -F * l / 2),
    },
  },
};

local fz_hinge_uz_rightward_1(F, l) = fz_common(F, l) + {
  name: 'fz_hinge_uz_rightward_1',

  dirichlet: qz_hinge_uz_rightward(0, l).dirichlet,
  elements: qz_hinge_uz_rightward(0, l).elements,

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(0) + test.My(3 * F * l),
      D: test.Fz(2 * F),
    },
    interpolation: {
      AB: test.Constant('My', 3 * F * l),
      BC: test.Linear('My', 3 * F * l, 2 * F * l),
      CD: test.Linear('My', 2 * F * l, 0),
    },
  },
};

local fz_hinge_uz_leftward_1(F, l) = fz_hinge_uz_rightward_1(F, l) + {
  name: 'fz_hinge_uz_leftward_1',
  elements: qz_hinge_uz_leftward(0, l).elements,

  expected: super.expected {
    interpolation: {
      BA: test.Constant('My', -3 * F * l),
      CB: test.Linear('My', -2 * F * l, -3 * F * l),
      DC: test.Linear('My', 0, -2 * F * l),
    },
  },
};

local fz_hinge_uz_rightward_2(F, l) = fz_hinge_uz_rightward_1(F, l) + {
  name: 'fz_hinge_uz_rightward_2',

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Uz'] }),
    CD: bvp.Frame2d(),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(F) + test.My(0),
      D: test.Fz(F),
    },
    interpolation: {
      AB: test.Linear('My', 0, F * l),
      BC: test.Constant('My', F * l),
      CD: test.Linear('My', F * l, 0),
    },
  },
};

local fz_hinge_uz_leftward_2(F, l) = fz_hinge_uz_rightward_2(F, l) + {
  name: 'fz_hinge_uz_leftward_2',
  elements: {
    BA: bvp.Frame2d(),
    CB: bvp.Frame2d(hinges={ B: ['Uz'] }),
    DC: bvp.Frame2d(),
  },

  expected: super.expected + {
    interpolation: {
      BA: test.Linear('My', -F * l, 0),
      CB: test.Constant('My', -F * l),
      DC: test.Linear('My', 0, -F * l),
    },
  },
};

local fz_hinge_phi_phi(F, l) = fz_common(F, l) + {
  name: 'fz_hinge_phi_phi',
  dirichlet: qz_hinge_phi_phi(0, l).dirichlet,
  elements: qz_hinge_phi_phi(0, l).elements,

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(F) + test.My(-F * l),
      D: test.Fz(F) + test.My(F * l),
    },
    interpolation: {
      AB: test.Linear('My', -F * l, 0),
      BC: test.Constant('My', 0),
      CD: test.Linear('My', 0, -F * l),
    },
  },
};

local fz_hinge_uz_phi_rightward_1(F, l) = fz_common(F, l) + {
  name: 'fz_hinge_uz_phi_rightward_1',
  dirichlet: qz_hinge_uz_phi_rightward(0, l).dirichlet,
  elements: qz_hinge_uz_phi_rightward(0, l).elements,

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(F) + test.My(-F * l),
      D: test.Fz(F) + test.My(F * l),
    },
    interpolation: {
      AB: test.Linear('My', -F * l, 0),
      BC: test.Constant('My', 0),
      CD: test.Linear('My', 0, -F * l),
    },
  },
};

local fz_hinge_uz_phi_leftward_1(F, l) = fz_hinge_uz_phi_rightward_1(F, l) + {
  name: 'fz_hinge_uz_phi_leftward_1',
  elements: qz_hinge_uz_phi_leftward(0, l).elements,

  expected: super.expected + {
    interpolation: {
      BA: test.Linear('My', 0, F * l),
      CB: test.Constant('My', 0),
      DC: test.Linear('My', F * l, 0),
    },
  },
};

local fz_hinge_uz_phi_rightward_2(F, l) = fz_common(F, l) + {
  name: 'fz_hinge_uz_phi_rightward_2',
  dirichlet: qz_hinge_uz_phi_rightward(0, l).dirichlet,

  elements: {
    AB: bvp.Frame2d(hinges={ B: ['Uz'] }),
    BC: bvp.Frame2d(hinges={ C: ['Phiy'] }),
    CD: bvp.Frame2d(),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(0) + test.My(F * l),
      D: test.Fz(2 * F) + test.My(2 * F * l),
    },
    interpolation: {
      AB: test.Constant('My', F * l),
      BC: test.Linear('My', F * l, 0),
      CD: test.Linear('My', 0, -2 * F * l),
    },
  },
};

local fz_hinge_uz_phi_leftward_2(F, l) = fz_hinge_uz_phi_rightward_2(F, l) + {
  name: 'fz_hinge_uz_phi_leftward_2',

  elements: {
    BA: bvp.Frame2d(hinges={ B: ['Uz'] }),
    CB: bvp.Frame2d(hinges={ C: ['Phiy'] }),
    DC: bvp.Frame2d(),
  },

  expected: super.expected + {
    interpolation: {
      BA: test.Constant('My', -F * l),
      CB: test.Linear('My', 0, -F * l),
      DC: test.Linear('My', 2 * F * l, 0),
    },
  },
};

local hinge_fails = {
  // We constrain most dofs, so that the provoked failure is always attributed to the hinge pattern.
  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    B: bvp.Uz(),
    C: bvp.Uz(),
    D: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
  },

  expected: {
    failure: 'unsupported.*hinge',
  },
};

local fz_hinge_uz_uz_fails(F, l) = fz_common(F, l) + hinge_fails {
  name: 'fz_hinge_uz_uz_fails',

  elements: {
    BA: bvp.Frame2d(),
    CB: bvp.Frame2d(hinges={ B: ['Uz'], C: ['Uz'] }),
    DC: bvp.Frame2d(),
  },
};

local fz_hinge_uz_phi_left_fails(F, l) = fz_common(F, l) + hinge_fails {
  name: 'fz_hinge_uz_phi_left_fails',

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Uz', 'Phiy'] }),
    CD: bvp.Frame2d(),
  },
};

local fz_hinge_uz_phi_right_fails(F, l) = fz_common(F, l) + hinge_fails {
  name: 'fz_hinge_uz_phi_right_fails',

  elements: {
    BA: bvp.Frame2d(),
    CB: bvp.Frame2d(hinges={ B: ['Uz', 'Phiy'] }),
    DC: bvp.Frame2d(),
  },
};

local fz_hinge_ux_ux_fails(F, l) = fz_common(F, l) + hinge_fails + {
  name: 'fz_hinge_ux_ux_fails',

  elements: {
    AB: bvp.Frame2d(),
    BC: bvp.Frame2d(hinges={ B: ['Ux'], C: ['Ux'] }),
    CD: bvp.Frame2d(),
  },
};

[
  qz_hinge_phi_rightward(1.5e3, 3),
  qz_hinge_phi_leftward(2.5e3, 4),
  qz_hinge_uz_rightward(0.5e3, 1.5),
  qz_hinge_uz_leftward(1.5e3, 0.5),
  qz_hinge_phi_phi(1.25e3, 2.75),
  qz_hinge_phi_phi_variation_1(1e3, 3.0),
  qz_hinge_phi_phi_variation_2(-2e3, 2.0),
  qz_hinge_phi_phi_variation_3(1.234e3, 2.0),
  qz_hinge_uz_phi_rightward(0.75e3, 1.5),
  qz_hinge_uz_phi_leftward(0.75e3, 1.5),
]
+
[
  fz_hinge_phi_rightward(12e3, 2.5),
  fz_hinge_phi_leftward(10e3, 2),
  fz_hinge_uz_rightward_1(10e3, 2),
  fz_hinge_uz_leftward_1(10e3, 2),
  fz_hinge_uz_rightward_2(10e3, 2),
  fz_hinge_uz_leftward_2(10e3, 2),
  fz_hinge_phi_phi(15e3, 2.5),
  fz_hinge_uz_phi_rightward_1(12.5e3, 1.5),
  fz_hinge_uz_phi_leftward_1(12.5e3, 1.5),
  fz_hinge_uz_phi_rightward_2(12.5e3, 1.75),
  fz_hinge_uz_phi_leftward_2(12.5e3, 1.75),
]
+
[
  fz_hinge_uz_uz_fails(10e3, 1),
  fz_hinge_uz_phi_left_fails(10e3, 1),
  fz_hinge_uz_phi_right_fails(10e3, 1),
  fz_hinge_ux_ux_fails(10e3, 1),
]
