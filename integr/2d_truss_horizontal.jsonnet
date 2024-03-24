local lib = import 'common.libsonnet';

local l = 1.0;
local default = lib.Defaults();

local single = {
  nodes: {
    A: [0, 0, 0],
    Z: [l, 0, 0],
  },

  elements: {
    AZ: default.Truss2d(),
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    Z: lib.Uz(),
  },
};

local multi = single {
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
    [elmt]: default.Truss2d()
    for elmt in ['AB', 'BC', 'CD', 'DE', 'EF', 'FZ']
  },

  dirichlet: super.dirichlet + {
    B: lib.Uz(),
    C: lib.Uz(),
    D: lib.Uz(),
    E: lib.Uz(),
    F: lib.Uz(),
  },
};

local nodal_fx = {
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

local single_nodal_fx = single + nodal_fx {
  name: 'single',
  description: 'Horizontal line, single truss',

  neumann: {
    Z: lib.Fx(100e3)
       // The vertical load doesn't cause any deformation, it goes straight into the support.
       // Specifying it here is to assert that the reaction force ends up in the results.
       + lib.Fz(-250e3),
  },

  expected: super.expected + {
    reaction: super.reaction + {
      Z: lib.Fz(250e3),
    },
  },
};

local single_element_fx = single {
  name: 'single_element_fx',
  description: 'Single horizontal truss, concentrated force',
  local fx = 123e3,
  local a = 0.7655 * l,

  neumann: {
    AZ: lib.Fx(-fx, a),
  },

  expected: {
    primary: {
      local eps = fx / (30000e6 * 0.01),
      local deltaU = eps * a,
      Z: lib.Ux(-deltaU),
    },
    reaction: {
      A: lib.Fx(fx),
    },
  },
};

local single_const_qx = single {
  name: 'single_const_qx',
  description: 'Single horizontal truss, const qx',
  local q = 10e3,

  neumann: {
    AZ: lib.qx(q),
  },

  expected: {
    reaction: {
      A: lib.Fx(-q * l),
    },
  },
};

local single_linear_qx = single {
  name: 'single_linear_qx',
  description: 'Single horizontal truss, linear qx',
  local q0 = -2.5e3,
  local q1 = 4.123e3,

  neumann: {
    AZ: lib.qx([q0, q1]),
  },

  expected: {
    reaction: {
      A: lib.Fx(-(q0 + q1) * l / 2),
    },
  },
};

local multi_const_qx = multi {
  name: 'multi_const_qx',
  description: 'Horizontal line, 6 trusses, constant element load',
  local q = 10e3,

  neumann: {
    [elmt]: lib.qx(q)
    for elmt in ['AB', 'BC', 'CD', 'DE', 'EF', 'FZ']
  },

  expected: {
    reaction: {
      A: lib.Fx(-q * l),
    },
  },
};

local multi_linear_qx = multi {
  name: 'multi_linear_qx',
  description: 'Horizontal line, 6 trusses, linear element load',
  local qe = 60e3,

  neumann: {
    // The load steps come from the x-coordinates of the nodes, they are not evenly spaced.
    AB: lib.qx([0 / 1 * qe, 1 / 5 * qe]),
    BC: lib.qx([1 / 5 * qe, 1 / 3 * qe]),
    CD: lib.qx([1 / 3 * qe, 1 / 2 * qe]),
    DE: lib.qx([1 / 2 * qe, 3 / 5 * qe]),
    EF: lib.qx([3 / 5 * qe, 3 / 4 * qe]),
    FZ: lib.qx([3 / 4 * qe, 1 / 1 * qe]),
  },

  expected: {
    reaction: {
      A: lib.Fx(-qe * l / 2),
    },
  },
};

local multi_fx = multi + nodal_fx {
  name: 'multi_fx',
  description: 'Horizontal line, 6 trusses, nodal force',
};

local single_ux = single {
  name: 'single_ux',
  description: 'Horizontal line, single truss ux prescribed',

  elements: {
    AZ: default.Truss2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    Z: lib.Ux(1 / 3 * 1e-3) + lib.Uz(),
  },

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

[
  single_nodal_fx,
  single_element_fx,
  single_const_qx,
  single_linear_qx,
  single_ux,
  multi_fx,
  multi_const_qx,
  multi_linear_qx,
  multi_ux,
]
