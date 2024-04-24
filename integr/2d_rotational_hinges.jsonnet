local lib = import 'common.libsonnet';

local default = lib.Defaults();

local tilted_l = {
  name: 'tilted_l',
  description: 'Example structure with nodal/element loads and rotational hinges',

  nodes: {
    A: [0, 0, 0],
    B: [3, 0, 0],
    C: [6, 0, 0],
    D: [6, 0, 3],
  },

  material: lib.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=0.01, Iyy=10e-6),

  neumann: {
    B: lib.Fz(3e3),
    BC: lib.Fz(4e3, 1.5),
    CD: lib.qz(-4e3),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz() + lib.Phiy(),
    D: lib.Ux() + lib.Uz(),
  },

  expected: {
    reaction: {
      A: lib.Fx(6e3) + lib.Fz(-1e3) + lib.My(3e3),
      D: lib.Fx(6e3) + lib.Fz(2e3),
    },
    interpolation: {
      AB: lib.Constant('Nx', -6e3) +
          lib.Constant('Vz', -1e3) +
          lib.Linear('My', 3e3, 0),
      BC: lib.Constant('Nx', -6e3) +
          lib.Constant('Vz', 2e3, range=[0, 1.5]) +
          lib.Constant('Vz', -2e3, range=[1.5, 3]) +
          lib.Linear('My', 0, 3e3, range=[0, 1.5]) +
          lib.Linear('My', 3e3, 0, range=[1.5, 3]),
      CD: lib.Constant('Nx', 2e3) +
          lib.Linear('Vz', -6e3, 6e3) +
          lib.Quadratic('My', eval=[[0, 0], [1.5, -4.5e3], [3, 0]]),
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

  material: lib.LinElast('default', E=210000e6, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=0.01, Iyy=10e-6),

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
    A: lib.Uz(),
    B: lib.Ux() + lib.Uz(),
    G: lib.Uz(),
    H: lib.Uz(),
  },

  neumann: {
    F: lib.Fx(2e3),
    BC: lib.Fz(4e3 * std.sin(30 * lib.pi / 180), 1.5) +
        lib.Fx(4e3 * std.cos(30 * lib.pi / 180), 1.5),
    CD: lib.qz(2e3),
  },

  expected: {
    tolerance: {
      polynomial: 5e-7,
    },
    // The sum of all horizontal forces:
    local Fh = 4e3 * std.cos(30 / 180.0 * lib.pi) + 2e3 + 2e3 * 3,

    reaction: {
      A: lib.Fz(1e3),
      B: lib.Fx(-Fh) + lib.Fz(-1e3),
      G: lib.Fz(15e3),
      H: lib.Fz(-7e3),
    },
    interpolation: {
      AB: lib.Constant('Nx', 0) +
          lib.Constant('Vz', 1e3) +
          lib.Linear('My', 0, 3e3),
      BC: lib.Constant('Nx', Fh, range=[0, 1.5]) +
          lib.Constant('Nx', 8e3, range=[1.5, 3]) +
          lib.Constant('Vz', 0, range=[0, 1.5]) +
          lib.Constant('Vz', -2e3, range=[1.5, 3]) +
          lib.Constant('My', 3e3, range=[0, 1.5]) +
          lib.Linear('My', 3e3, 0, range=[1.5, 3]),
      local lcd = 3 * std.sqrt(2),
      CD: lib.Constant('Nx', 5e3 * std.sqrt(2)) +
          lib.Linear('Vz', 2e3 * lcd / 2, -2e3 * lcd / 2) +
          lib.Quadratic('My', eval=[[0, 0], [0.5 * lcd, 2e3 * lcd * lcd / 8], [lcd, 0]]),
      DE: lib.Constant('Nx', 2e3) +
          lib.Constant('Vz', -8e3) +
          lib.Linear('My', 0, -12e3),
      EF: lib.Constant('Nx', 0) +
          lib.Constant('Vz', 2e3) +
          lib.Linear('My', -3e3, 0),
      EG: lib.Constant('Nx', 0) +
          lib.Constant('Vz', -8e3) +
          lib.Linear('My', -9e3, -21e3),
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
