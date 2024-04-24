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
    interpolation: {
      AZ: lib.Constant('Nx', 100e3) +
          lib.Linear('Ux', 0, 1 / 3 * 1e-3),
    },
  },
};

local single_element_fx = single {
  name: 'single_element_fx',
  description: 'Single horizontal truss, concentrated force',
  local fx = -123e3,
  local a = 0.7655 * l,

  neumann: {
    AZ: lib.Fx(fx, a),
  },

  expected: {
    local eps = fx / (30000e6 * 0.01),
    local deltaU = eps * a,

    primary: {
      Z: lib.Ux(deltaU),
    },
    reaction: {
      A: lib.Fx(-fx),
    },
    interpolation: {
      AZ: lib.Constant('Nx', fx, range=[0, a]) +
          lib.Constant('Nx', 0, range=[a, l]) +
          lib.Linear('Ux', 0, deltaU, range=[0, a]) +
          lib.Constant('Ux', deltaU, range=[a, l]),
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
    interpolation: {
      AZ: lib.Linear('Nx', q * l, 0) +
          lib.Quadratic('Ux', eval=[[0, 0], [l, 0.5 * q * l * l / (30000e6 * 0.01)]]),
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
    local H = (q0 + q1) * l / 2,

    reaction: {
      A: lib.Fx(-H),
    },
    interpolation: {
      AZ: lib.Quadratic('Nx', eval=[[0, H], [l, 0]]) +
          lib.Cubic('Ux', eval=[[0, 0]]),
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
    interpolation: {
      local ql = q * l,
      AB: lib.Linear('Nx', ql, 4 / 5 * ql),
      BC: lib.Linear('Nx', 4 / 5 * ql, 2 / 3 * ql),
      CD: lib.Linear('Nx', 2 / 3 * ql, ql / 2),
      DE: lib.Linear('Nx', ql / 2, 2 / 5 * ql),
      EF: lib.Linear('Nx', 2 / 5 * ql, 1 / 4 * ql),
      FZ: lib.Linear('Nx', 1 / 4 * ql, 0),
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
    interpolation: {
      local N(x) = qe * l / 2 - qe / (2 * l) * x * x,
      local samples(x0, x1, n) = std.map(function(eval) [eval[0] - x0, eval[1]], lib.Samples(N, x0, x1, n)),
      AB: lib.Quadratic('Nx', eval=samples(0 / 1 * l, 1 / 5 * l, 5)),
      BC: lib.Quadratic('Nx', eval=samples(1 / 5 * l, 1 / 3 * l, 5)),
      CD: lib.Quadratic('Nx', eval=samples(1 / 3 * l, 1 / 2 * l, 5)),
      DE: lib.Quadratic('Nx', eval=samples(1 / 2 * l, 3 / 5 * l, 5)),
      EF: lib.Quadratic('Nx', eval=samples(3 / 5 * l, 3 / 4 * l, 5)),
      FZ: lib.Quadratic('Nx', eval=samples(3 / 4 * l, 1 / 1 * l, 5)),
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

  local deltaU = 1 / 3 * 1e-3,

  dirichlet: super.dirichlet + {
    Z: lib.Ux(deltaU) + lib.Uz(),
  },

  neumann: {},

  expected: {
    interpolation: {
      local EA = 0.01 * 30000e6,
      local eps = deltaU / l,
      [elmt]: lib.Constant('Nx', eps * EA)
      for elmt in ['AB', 'BC', 'CD', 'DE', 'EF', 'FZ']
    },
  },
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
