local lib = import 'common.libsonnet';

local default = lib.Defaults();

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

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

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
    A: lib.Ux() + lib.Uz(),
    E: lib.Uz(),
  },

  neumann: {
    F: lib.Fz(-5e3),
    G: lib.Fz(-10e3),
    H: lib.Fz(-10e3),
    I: lib.Fz(-10e3),
    J: lib.Fz(-5e3),
  },

  expected: {
    reaction: {
      A: lib.Fx(0) + lib.Fz(20e3),
      E: lib.Fz(20e3),
    },
  },
};

[bridge(1), bridge(10), bridge(25)]
