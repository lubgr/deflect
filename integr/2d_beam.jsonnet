local lib = import 'common.libsonnet';

local default = lib.Defaults();

local common(l, E, Iyy) = {
  material: lib.LinElast('default', E=E, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=1, Iyy=Iyy),

  elements: {
    AB: default.Frame2d(),
  },
};

local horizontal(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Uz(),
  },
};

local vertical(l) = {
  nodes: {
    A: [0, 0, 0],
    B: [0, 0, l],
  },

  dirichlet: {
    A: lib.Ux() + lib.Uz(),
    B: lib.Ux(),
  },
};

local with_const_qz(q, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) {
  name: 'element_const_qz_%g' % q,

  neumann: {
    AB: lib.qz(q),
  },

  expected: {
    local EI = E * Iyy,
    local l3 = std.pow(l, 3),
    local phiy = q * l3 / (24 * EI),

    primary: {
      A: lib.Phiy(phiy),
      B: lib.Phiy(-phiy),
    },
    reaction: {
      A: lib.Fz(q * l / 2),
      B: lib.Fz(q * l / 2),
    },
    interpolation: {
      AB: lib.Constant('Nx', 0) +
          lib.Linear('Vz', q * l / 2, -q * l / 2) +
          lib.Quadratic('My', eval=[[0, 0], [l / 2, q * l * l / 8], [l, 0]]) +
          lib.Cubic('Phiy', eval=[[0, -phiy], [l / 2, 0], [l, phiy]]) +
          lib.Quartic('Uz', eval=lib.Samples(
            function(x)
              local l4 = std.pow(l, 4);
              q * l4 / (24 * EI) * (x / l - 2 * std.pow(x, 3) / l3 + std.pow(x, 4) / l4), 0, l, 20
          )),
    },
  },
};

local with_const_qz_vertical(q, l, E, Iyy) = with_const_qz(q, l, E, Iyy) + vertical(l) + {
  name: 'element_const_qz_%g_vertical' % q,

  expected: super.expected + {
    primary: super.primary,
    reaction: {
      A: lib.Fx(-q * l / 2),
      B: lib.Fx(-q * l / 2),
    },
    interpolation: super.interpolation,
  },
};

local with_linear_qz(q0, q1, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'element_linear_qz_%.1f_%.1f' % [q0, q1],

  neumann: {
    AB: lib.qz([q0, q1]),
  },

  expected: {
    primary: {
      A: lib.Phiy((8 * q0 + 7 * q1) * std.pow(l, 3) / (360 * E * Iyy)),
      B: lib.Phiy(-(7 * q0 + 8 * q1) * std.pow(l, 3) / (360 * E * Iyy)),
    },
    local Av = (2 * q0 + q1) * l / 6,
    local Bv = (q0 + 2 * q1) * l / 6,
    reaction: {
      A: lib.Fz(Av),
      B: lib.Fz(Bv),
    },
    interpolation: {
      local my(x) =
        local aux(xbar) = xbar / l - std.pow(xbar / l, 3);
        l * l / 6 * (q0 * aux(l - x) + q1 * aux(x)),
      local uz(x) =
        local EI = E * Iyy;
        local aux(xbar) = (xbar / l - std.pow(xbar / l, 3)) * (7 - 3 * std.pow(xbar / l, 2));
        std.pow(l, 4) / (360 * EI) * (aux(l - x) * q0 + aux(x) * q1),

      AB: lib.Constant('Nx', 0) +
          (
            if q0 == q1 then
              lib.Linear('Vz', Av, -Bv) +
              lib.Quadratic('My', eval=lib.Samples(my, 0, l, 10)) +
              lib.Quartic('Uz', eval=lib.Samples(uz, 0, l, 10))
            else
              lib.Quadratic('Vz', eval=[[0, Av], [l, -Bv]]) +
              lib.Cubic('My', eval=lib.Samples(my, 0, l, 10)) +
              lib.Quintic('Uz', eval=lib.Samples(uz, 0, l, 10))
          ),
    },
  },
};

local with_element_fz(F, a, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'element_fz_%.1f_%.2f' % [F, a],

  neumann: {
    AB: lib.Fz(F, a),
  },

  expected: {
    primary: {
      A: lib.Phiy(F * a * (l - a) * (2 * l - a) / (6 * E * Iyy * l)),
      B: lib.Phiy(-F * a * (l - a) * (l + a) / (6 * E * Iyy * l)),
    },

    local Av = F * (1 - a / l),
    local Bv = F * a / l,

    reaction: {
      A: lib.Fz(Av),
      B: lib.Fz(Bv),
    },

    interpolation: {
      local b = l - a,
      local EI = E * Iyy,
      local leftUz(x) = F * a * b * b / (6 * EI) * (
        (1 + l / b) * x / l - x * x * x / (a * b * l)
      ),
      local rightUz(x) = F * a * a * b / (6 * EI) * (
        (1 + l / a) * (l - x) / l - std.pow(l - x, 3) / (a * b * l)
      ),

      AB: lib.Constant('Nx', 0) +
          (
            if a == 0 || a == l then
              lib.Constant('Vz', 0) +
              lib.Constant('My', 0) +
              lib.Constant('Uz', 0)
            else
              lib.Constant('Vz', Av, range=[0, a]) +
              lib.Constant('Vz', -Bv, range=[a, l]) +
              lib.Linear('My', 0, F * a * (1 - a / l), range=[0, a]) +
              lib.Linear('My', F * a * (1 - a / l), 0, range=[a, l]) +
              lib.Cubic('Uz', range=[0, a], eval=lib.Samples(leftUz, 0, a, 5)) +
              lib.Cubic('Uz', range=[a, l], eval=lib.Samples(rightUz, a, l, 5))
          ),
    },
  },
};

local with_element_fz_vertical(F, a, l, E, Iyy) = with_element_fz(F, a, l, E, Iyy) + vertical(l) + {
  name: 'element_fz_%g_%.2f_vertical' % [F, a],

  expected: super.expected + {
    primary: super.primary,
    reaction: {
      A: lib.Fx(-F * (1 - a / l)),
      B: lib.Fx(-F * a / l),
    },
    interpolation: super.interpolation,
  },
};

local with_element_my(M, a, l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'element_my_%.1f_%.2f' % [M, a],

  neumann: {
    AB: lib.My(M, a),
  },

  expected: {
    local alpha = a / l,
    local beta = (l - a) / l,
    local EI = E * Iyy,
    local leftPhiy = M * (1 - 3 * beta * beta) / (6 * EI / l),
    local rightPhiy = M * (1 - 3 * alpha * alpha) / (6 * EI / l),

    primary: {
      A: lib.Phiy(leftPhiy),
      B: lib.Phiy(rightPhiy),
    },
    reaction: {
      A: lib.Fz(M / l),
      B: lib.Fz(-M / l),
    },
    interpolation: {
      local leftUz(x) =
        local xi = x / l;
        xi * M * l * l / (6 * EI) * (1 - 3 * beta * beta - xi * xi),
      local rightUz(x) =
        local xib = (l - x) / l;
        -xib * M * l * l / (6 * EI) * (1 - 3 * alpha * alpha - xib * xib),

      AB: lib.Constant('Nx', 0) +
          lib.Constant('Vz', M / l) +
          (
            if a == 0 then
              lib.Linear('My', -M, 0) +
              lib.Quadratic('Phiy', eval=[[0, -leftPhiy], [l, -rightPhiy]]) +
              lib.Cubic('Uz', eval=lib.Samples(rightUz, 0, l, 10))
            else if a == l then
              lib.Linear('My', 0, M) +
              lib.Quadratic('Phiy', eval=[[0, -leftPhiy], [l, -rightPhiy]]) +
              lib.Cubic('Uz', eval=lib.Samples(leftUz, 0, l, 10))
            else
              lib.Linear('My', 0, M * a / l, range=[0, a]) +
              lib.Linear('My', M * (a / l - 1), 0, range=[a, l]) +
              lib.Quadratic('Phiy', range=[0, a], eval=[[0, -leftPhiy]]) +
              lib.Quadratic('Phiy', range=[a, l], eval=[[l, -rightPhiy]]) +
              lib.Cubic('Uz', range=[0, a], eval=lib.Samples(leftUz, 0, a, 5)) +
              lib.Cubic('Uz', range=[a, l], eval=lib.Samples(rightUz, a, l, 5))
          ),
    },
  },
};

local with_all(l, E, Iyy) = common(l, E, Iyy) + horizontal(l) + {
  name: 'all_loads_01',

  local fconcentrated = with_element_fz(11.5e3, 0.3 * l, l, E, Iyy),
  local myconcentrated = with_element_my(2.5e3, 0.7 * l, l, E, Iyy),
  local qconst = with_const_qz(1.25e3, l, E, Iyy),
  local qlinear = with_linear_qz(-1.25e3, 1.25e3, l, E, Iyy),

  neumann: {
    AB: fconcentrated.neumann.AB + myconcentrated.neumann.AB + qconst.neumann.AB + qlinear.neumann.AB,
  },

  expected: {
    primary: {
      local phiy(from, node) = from.expected.primary[node][0].Phiy,
      A: lib.Ux(0) + lib.Uz(0) +
         lib.Phiy(phiy(fconcentrated, 'A') + phiy(myconcentrated, 'A') + phiy(qconst, 'A') + phiy(qlinear, 'A')),
      B: lib.Ux(0) + lib.Uz(0) +
         lib.Phiy(phiy(fconcentrated, 'B') + phiy(myconcentrated, 'B') + phiy(qconst, 'B') + phiy(qlinear, 'B')),
    },
    reaction: {
      local fz(from, node) = from.expected.reaction[node][0].Fz,
      A: lib.Fx(0) +
         lib.Fz(fz(fconcentrated, 'A') + fz(myconcentrated, 'A') + fz(qconst, 'A') + fz(qlinear, 'A')),
      B: lib.Fz(fz(fconcentrated, 'B') + fz(myconcentrated, 'B') + fz(qconst, 'B') + fz(qlinear, 'B')),
    },
  },
};

[
  with_const_qz(q=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_const_qz(q=2.5e3, l=6, E=12345e6, Iyy=8e-6),
  with_const_qz_vertical(q=2.5e3, l=6, E=12345e6, Iyy=8e-6),

  with_linear_qz(q0=1e3, q1=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=0e3, q1=1e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=1e3, q1=0e3, l=1, E=10000e6, Iyy=10e-6),
  with_linear_qz(q0=-2.5e3, q1=10e3, l=2, E=12345e6, Iyy=8e-6),
  with_linear_qz(q0=2.5e3, q1=-8e3, l=2, E=12345e6, Iyy=8e-6),
]
+
[
  with_element_fz(F=10e3, a=a, l=4, E=10000e6, Iyy=10e-6)
  for a in [0, 0.25, 2, 3.9, 4]
]
+
[
  with_element_fz_vertical(F=10e3, a=a, l=4, E=10000e6, Iyy=10e-6)
  for a in [0, 0.25, 2, 3.9, 4]
]
+
[
  with_element_my(M=2.5e3, a=a, l=10, E=10000e6, Iyy=10e-6)
  for a in [0, 10 / 3, 5, 10]
]
+
[
  with_all(l=4, E=10000e6, Iyy=10e-6),
]
