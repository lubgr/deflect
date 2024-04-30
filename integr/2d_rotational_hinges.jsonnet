local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

local tilted_l = {
  name: 'tilted_l',
  description: 'Example structure with nodal/element loads and rotational hinges',

  nodes: {
    A: [0, 0, 0],
    B: [3, 0, 0],
    C: [6, 0, 0],
    D: [6, 0, 3],
  },

  material: bvp.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=10e-6),

  neumann: {
    B: bvp.Fz(3e3),
    BC: bvp.Fz(4e3, 1.5),
    CD: bvp.qz(-4e3),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz() + bvp.Phiy(),
    D: bvp.Ux() + bvp.Uz(),
  },

  expected: {
    reaction: {
      A: test.Fx(6e3) + test.Fz(-1e3) + test.My(3e3),
      D: test.Fx(6e3) + test.Fz(2e3),
    },
    interpolation: {
      AB: test.Constant('Nx', -6e3) +
          test.Constant('Vz', -1e3) +
          test.Linear('My', 3e3, 0),
      BC: test.Constant('Nx', -6e3) +
          test.Constant('Vz', 2e3, range=[0, 1.5]) +
          test.Constant('Vz', -2e3, range=[1.5, 3]) +
          test.Linear('My', 0, 3e3, range=[0, 1.5]) +
          test.Linear('My', 3e3, 0, range=[1.5, 3]),
      CD: test.Constant('Nx', 2e3) +
          test.Linear('Vz', -6e3, 6e3) +
          test.Quadratic('My', eval=[[0, 0], [1.5, -4.5e3], [3, 0]]),
    },
  },
};

local tilted_l_hinges1 = tilted_l {
  elements: {
    AB: default.Frame2d(hinges={ B: ['Phiy'] }),
    BC: default.Frame2d(hinges={ C: ['Phiy'] }),
    CD: default.Frame2d(hinges={}),
  },
};

local tilted_l_hinges2 = tilted_l {
  elements: {
    AB: default.Frame2d(),
    BC: default.Frame2d(hinges={ B: ['Phiy'] }),
    CD: default.Frame2d(hinges={ C: ['Phiy'] }),
  },
};

local tilted_l_hinges3 = tilted_l {
  elements: {
    AB: default.Frame2d(),
    BC: default.Frame2d(hinges={ B: ['Phiy'], C: ['Phiy'] }),
    CD: default.Frame2d(),
  },
};

local double_plateau = {
  name: 'double_plateau',
  description: 'Example structure with nodal/element loads and rotational hinges',

  nodes: {
    A: [0, 0, 0],
    B: [3, 0, 0],
    C: [6, 0, 0],
    D: [9, 0, 3],
    E: [10.5, 0, 3],
    F: [10.5, 0, 4.5],
    G: [12, 0, 3],
    H: [15, 0, 3],
  },

  material: bvp.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=10e-6),

  elements: {
    AB: default.Frame2d(),
    BC: default.Frame2d(),
    CD: default.Frame2d(),
    DE: default.Frame2d(),
    EF: default.Frame2d(),
    EG: default.Frame2d(),
    GH: default.Frame2d(),
  },

  dirichlet: {
    A: bvp.Uz(),
    B: bvp.Ux() + bvp.Uz(),
    G: bvp.Uz(),
    H: bvp.Uz(),
  },

  neumann: {
    F: bvp.Fx(2e3),
    BC: bvp.Fz(4e3 * std.sin(30 * bvp.pi / 180), 1.5) +
        bvp.Fx(4e3 * std.cos(30 * bvp.pi / 180), 1.5),
    CD: bvp.qz(2e3),
  },

  expected: {
    tolerance: {
      polynomial: 5e-7,
    },
    // The sum of all horizontal forces:
    local Fh = 4e3 * std.cos(30 / 180.0 * bvp.pi) + 2e3 + 2e3 * 3,

    reaction: {
      A: test.Fz(1e3),
      B: test.Fx(-Fh) + test.Fz(-1e3),
      G: test.Fz(15e3),
      H: test.Fz(-7e3),
    },
    interpolation: {
      AB: test.Constant('Nx', 0) +
          test.Constant('Vz', 1e3) +
          test.Linear('My', 0, 3e3),
      BC: test.Constant('Nx', Fh, range=[0, 1.5]) +
          test.Constant('Nx', 8e3, range=[1.5, 3]) +
          test.Constant('Vz', 0, range=[0, 1.5]) +
          test.Constant('Vz', -2e3, range=[1.5, 3]) +
          test.Constant('My', 3e3, range=[0, 1.5]) +
          test.Linear('My', 3e3, 0, range=[1.5, 3]),
      local lcd = 3 * std.sqrt(2),
      CD: test.Constant('Nx', 5e3 * std.sqrt(2)) +
          test.Linear('Vz', 2e3 * lcd / 2, -2e3 * lcd / 2) +
          test.Quadratic('My', eval=[[0, 0], [0.5 * lcd, 2e3 * lcd * lcd / 8], [lcd, 0]]),
      DE: test.Constant('Nx', 2e3) +
          test.Constant('Vz', -8e3) +
          test.Linear('My', 0, -12e3),
      EF: test.Constant('Nx', 0) +
          test.Constant('Vz', 2e3) +
          test.Linear('My', -3e3, 0),
      EG: test.Constant('Nx', 0) +
          test.Constant('Vz', -8e3) +
          test.Linear('My', -9e3, -21e3),
    },
  },
};

local double_plateau_hinges1 = double_plateau {
  elements: super.elements + {
    CD: default.Frame2d(hinges={ C: ['Phiy'], D: ['Phiy'] }),
  },
};

local double_plateau_hinges2 = double_plateau {
  elements: super.elements + {
    BC: default.Frame2d(hinges={ C: ['Phiy'] }),
    CD: default.Frame2d(hinges={ D: ['Phiy'] }),
  },
};

local double_plateau_hinges3 = double_plateau {
  elements: super.elements + {
    BC: default.Frame2d(hinges={ C: ['Phiy'] }),
    DE: default.Frame2d(hinges={ D: ['Phiy'] }),
  },
};

[
  // We test all possible arrangement of the hinges:
  tilted_l_hinges1,
  tilted_l_hinges2,
  tilted_l_hinges3,
  double_plateau_hinges1,
  double_plateau_hinges2,
  double_plateau_hinges3,
]
