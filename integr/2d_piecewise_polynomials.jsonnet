local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local multi_load_interval_inclined_beam(q0, q1, F, M, l, angle) = {
  name: 'multi_load_interval_inclined_beam',
  description: 'Inclined beam with multiple element loads',

  local alpha = angle * bvp.pi / 180.0,

  nodes: {
    A: [0, 0, 0],
    B: [l * std.cos(alpha), 0, l * std.sin(alpha)],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: bvp.Frame2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
  },

  links: {
    B: bvp.InclinedSupportUxUz(alpha),
  },

  neumann: {
    AB: bvp.qz([q0, q1]) +
        bvp.Fz(F[0], 1 / 4 * l) +
        bvp.My(M, 5 / 12 * l) +
        bvp.Fz(F[1], 7 / 12 * l) +
        bvp.Fz(F[2], 11 / 12 * l),
  },

  expected: {
    interpolation: {
      local A = F[0] * 3 / 4 + F[1] * 5 / 12 + F[2] / 12 + M / l + q0 * l / 3 + q1 * l / 6,
      local vz0(x) = A - q0 * x - (q1 - q0) / (2 * l) * x * x,
      local vz1(x) = vz0(x) - F[0],
      local my0(x) = A * x - q0 / 2 * x * x - (q1 - q0) / (6 * l) * std.pow(x, 3),
      local my1(x) = my0(x) - F[0] * (x - l / 4),
      local my2(x) = my1(x) - M,
      local vz3(x) = vz1(x) - F[1],
      local my3(x) = my2(x) - F[1] * (x - 7 / 12 * l),
      local vz4(x) = vz3(x) - F[2],
      local my4(x) = my3(x) - F[2] * (x - 11 / 12 * l),

      AB: test.Quadratic('Vz', range=[0, l / 4], eval=test.Samples(vz0, 0, l / 4, 5)) +
          test.Quadratic('Vz', range=[l / 4, 7 / 12 * l], eval=test.Samples(vz1, l / 4, 5 / 12 * l, 5)) +
          test.Cubic('My', range=[0, l / 4], eval=test.Samples(my0, 0, l / 4, 5)) +
          test.Cubic('My', range=[l / 4, 5 / 12 * l], eval=test.Samples(my1, l / 4, 5 / 12 * l, 5)) +
          test.Cubic('My', range=[5 / 12 * l, 7 / 12 * l], eval=test.Samples(my2, 5 / 12 * l, 7 / 12 * l, 5)) +
          test.Quadratic('Vz', range=[7 / 12 * l, 11 / 12 * l], eval=test.Samples(vz3, 7 / 12 * l, 11 / 12 * l, 5)) +
          test.Cubic('My', range=[7 / 12 * l, 11 / 12 * l], eval=test.Samples(my3, 7 / 12 * l, 11 / 12 * l, 5)) +
          test.Quadratic('Vz', range=[11 / 12 * l, l], eval=test.Samples(vz4, 11 / 12 * l, l, 5)) +
          test.Cubic('My', range=[11 / 12 * l, l], eval=test.Samples(my4, 11 / 12 * l, l, 5)),
    },
  },
};

local multi_torque_beam = {
  name: 'multi_torque_beam',

  local q = 1e3,
  local M = 2.5e3,
  local l = 10,

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: bvp.Frame2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    B: bvp.Uz(),
  },

  neumann: {
    AB: bvp.qz(q) +
        bvp.My(M, 1 * l / 6) +
        bvp.My(M, 2 * l / 6) +
        bvp.My(M, 3 * l / 6) +
        bvp.My(M, 4 * l / 6) +
        bvp.My(M, 5 * l / 6),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(q * l / 2 + 5 * M / l),
      B: test.Fz(q * l / 2 - 5 * M / l),
    },
    interpolation: {
      // Ensures Vz is reduced to a single polynomial over the entire length
      AB: test.Linear('Vz', q * l / 2 + 5 * M / l, -q * l / 2 + 5 * M / l, range=[0, l]),
    },
  },
};

[
  multi_load_interval_inclined_beam(q0=1e3, q1=2e3, F=[2e3, 3e3, 4e3], M=0.5e3, l=5, angle=45),
  multi_load_interval_inclined_beam(q0=0, q1=1e3, F=[-1e3, -2e3, -3e3], M=-1e3, l=2, angle=20),
  multi_load_interval_inclined_beam(q0=-5e3, q1=2e3, F=[1e3, -2e3, 1e3], M=1e3, l=3.5, angle=123.456),
]
+
[
  multi_torque_beam,
]
