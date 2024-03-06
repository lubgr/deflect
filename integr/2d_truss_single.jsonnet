local lib = import 'common.libsonnet';

local l = 1.0;
local default = lib.Defaults();

local common = {
  nodes: {
    A: [0, 0, 0],
    Z: [l, 0, 0],
  },

  material: lib.LinElast('default', E=30000 * 1e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    Z: lib.Uz(),
  },

  neumann: {
    Z: lib.Fx(100e3),
  },

  expected: {
    primary: {
      A: lib.Ux() + lib.Uz(),
      Z: lib.Ux(1 / 3 * 1e-3) + lib.Uz(),
    },
    reaction: {
      A: lib.Fx(-100e3),
    },
  },
};

local single = common {
  name: 'single',
  description: 'Horizontal line, single truss',

  neumann: {
    Z: lib.Fx(100e3)
       // The vertical load doesn't cause any deformation, it goes straight into the support.
       // Specifying it here is to assert that the reaction force ends up in the results.
       + lib.Fz(-250e3),
  },

  elements: {
    AZ: default.Truss2d(),
  },

  expected: super.expected + {
    reaction: super.reaction + {
      Z: lib.Fz(250e3),
    },
  },
};

local multi = common {
  name: 'multi',
  description: 'Horizontal line, 6 trusses',

  // Semi-random nodes to partition the line into several sub-elements. The resulting assertions can
  // stay identical between then single-element and the multi-element line.
  nodes: super.nodes + {
    B: [l / 5, 0, 0],
    C: [l / 3, 0, 0],
    D: [l / 2, 0, 0],
    E: [3 * l / 5, 0, 0],
    F: [3 * l / 4, 0, 0],
  },

  elements: {
    AB: default.Truss2d(),
    BC: default.Truss2d(),
    CD: default.Truss2d(),
    DE: default.Truss2d(),
    EF: default.Truss2d(),
    FZ: default.Truss2d(),
  },

  dirichlet: super.dirichlet + {
    B: lib.Uz(),
    C: lib.Uz(),
    D: lib.Uz(),
    E: lib.Uz(),
    F: lib.Uz(),
  },
};

local single_ux = single {
  name: 'single_ux',
  description: 'Horizontal line, single truss ux prescribed',

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    Z: lib.Ux(1 / 3 * 1e-3) + lib.Uz(),
  },

  neumann: {},

  expected: {
    failure: '.*all.*degrees.*of.*freedom.*Dirichlet.*',
  },
};

local multi_ux = multi {
  name: 'multi_ux',
  description: 'Horizontal line, 6 trusses ux prescribed',

  dirichlet: super.dirichlet + {
    Z: lib.Ux(1 / 3 * 1e-3) + lib.Uz(),
  },

  neumann: {},
};

[single, multi, single_ux, multi_ux]
