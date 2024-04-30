local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

local bridge(a) = {
  name: 'bridge_%.2f' % a,
  description: 'Truss bridge with vertical nodal loads',

  nodes: {
    A: [0, 0, 0],
    B: [a, 0, 0],
    C: [2 * a, 0, 0],
    D: [3 * a, 0, 0],
    E: [4 * a, 0, 0],

    F: [0, 0, a],
    G: [a, 0, a],
    H: [2 * a, 0, a],
    I: [3 * a, 0, a],
    J: [4 * a, 0, a],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  local trusses(ids) = {
    [id]: default.Truss2d()
    for id in ids
  },

  elements: trusses([
    'AB',
    'BC',
    'CD',
    'DE',
    'FG',
    'GH',
    'HI',
    'IJ',
    'AF',
    'BG',
    'CH',
    'DI',
    'EJ',
    'FB',
    'GC',
    'CI',
    'DJ',
  ]),

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    E: bvp.Uz(),
  },

  neumann: {
    F: bvp.Fz(-5e3),
    G: bvp.Fz(-10e3),
    H: bvp.Fz(-10e3),
    I: bvp.Fz(-10e3),
    J: bvp.Fz(-5e3),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(20e3),
      E: test.Fz(20e3),
    },
    interpolation: {
      local Nx(value) = test.Constant('Nx', value),
      // Lower horizontal:
      AB: Nx(0),
      BC: Nx(15e3),
      CD: Nx(15e3),
      DE: Nx(0),
      // Upper horizontal:
      FG: Nx(-15e3),
      GH: Nx(-20e3),
      HI: Nx(-20e3),
      IJ: Nx(-15e3),
      // Vertical:
      AF: Nx(-20e3),
      BG: Nx(-15e3),
      CH: Nx(-10e3),
      DI: Nx(-15e3),
      EJ: Nx(-20e3),
      // Diagonal:
      FB: Nx(15e3 * std.sqrt(2)),
      GC: Nx(5e3 * std.sqrt(2)),
      CI: Nx(5e3 * std.sqrt(2)),
      DJ: Nx(15e3 * std.sqrt(2)),
    },
  },
};

[bridge(1), bridge(10), bridge(25)]
