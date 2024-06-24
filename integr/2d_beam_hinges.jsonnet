local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local beams_with_shear_hinge(q, l, E, Iyy, hinge) = {
  name: 'beams_shear_hinge',

  nodes: {
    A: [0, 0, 0],
    B: [2 * l, 0, 0],
    C: [3 * l, 0, 0],
  },

  material: bvp.LinElast('default', E=E, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=Iyy, Izz=10e-6),

  elements: {
    assert hinge == 'left' || hinge == 'right',
    AB: bvp.Frame2d(hinges=if hinge == 'left' then { B: ['Uz'] } else {}),
    BC: bvp.Frame2d(hinges=if hinge == 'right' then { B: ['Uz'] } else {}),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    C: bvp.Uz() + bvp.Phiy(),
  },

  neumann: {
    AB: bvp.qz(q),
    BC: bvp.qz(q),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(2 * q * l),
      C: test.Fz(q * l) + test.My(-3 / 2 * q * l * l),
    },
    interpolation: {
      local EI = E * Iyy,
      local l3 = std.pow(l, 3),
      local phiyAB(x) = 1 / EI * (q * l * x * x - q / 6 * std.pow(x, 3) - 9 / 2 * q * l3),
      local uzAB(x) = -1 / EI * (q * l / 3 * std.pow(x, 3) - q / 24 * std.pow(x, 4)
                                 - 9 / 2 * q * l3 * x),
      local phiyBC(x) = 1 / EI * (2 * q * l * l * x - q / 6 * std.pow(x, 3) - 11 / 6 * q * l3),
      local uzBC(x) = -1 / EI * (q * std.pow(l * x, 2) - q / 24 * std.pow(x, 4)
                                 - 11 / 6 * q * l3 * x + 7 / 8 * q * std.pow(l, 4)),

      // The hinge decouples the two degrees of freedom here. Depending on the
      // parameterisation of length, load and material properties, there could be
      // scenarios where this assertion fails, but that's very very unlikely.
      assert uzAB(2 * l) != uzBC(0),

      AB: test.Linear('Vz', 2 * q * l, 0) +
          test.Quadratic('My', eval=[[0, 0], [2 * l, 2 * q * l * l]]) +
          test.Cubic('Phiy', eval=test.Samples(phiyAB, 0, 2 * l, 7)) +
          test.Quartic('Uz', eval=test.Samples(uzAB, 0, 2 * l, 7)),
      BC: test.Linear('Vz', 0, -q * l) +
          test.Quadratic('My', eval=[[0, 2 * q * l * l], [l, 3 / 2 * q * l * l]]) +
          test.Cubic('Phiy', eval=test.Samples(phiyBC, 0, l, 5)) +
          // Ensures that the interpolation decouples the nodal displacement from the
          // one computed on the element level because of the hinge.
          test.Quartic('Uz', eval=test.Samples(uzBC, 0, l, 5)),
    },
  },
};

local beams_with_rot_hinge(q, l, E, Iyy, hinge) = {
  name: 'beams_rot_hinge',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [2 * l, 0, 0],
  },

  material: bvp.LinElast('default', E=E, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=0.01, Iyy=Iyy, Izz=10e-6),

  elements: {
    assert hinge == 'left' || hinge == 'right',
    AB: bvp.Frame2d(hinges=if hinge == 'left' then { B: ['Phiy'] } else {}),
    BC: bvp.Frame2d(hinges=if hinge == 'right' then { B: ['Phiy'] } else {}),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Uz(),
    C: bvp.Uz() + bvp.Phiy(),
  },

  neumann: {
    AB: bvp.qz(q),
    BC: bvp.qz(q),
  },

  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(q * l / 2),
      C: test.Fz(3 / 2 * q * l) + test.My(q * l * l),
    },
    interpolation: {
      local EI = E * Iyy,
      local l3 = std.pow(l, 3),
      local l4 = std.pow(l, 4),
      local phiyAB(x) = 1 / EI * (q * l / 4 * x * x - q / 6 * std.pow(x, 3) - q * l3 / 3),
      local uzAB(x) = -1 / EI * (q * l / 12 * std.pow(x, 3) - q / 24 * std.pow(x, 4)
                                 - q * l3 / 3 * x),
      local phiyBC(x) = 1 / EI * (-q * l / 4 * x * x - q / 6 * std.pow(x, 3) + 5 / 12 * q * l3),
      local uzBC(x) = -1 / EI * (-q * l / 12 * std.pow(x, 3) - q / 24 * std.pow(x, 4)
                                 + 5 / 12 * q * l3 * x - 7 / 24 * q * l4),

      AB: test.Linear('Vz', q * l / 2, -q * l / 2) +
          test.Quadratic('My', eval=[[0, 0], [l, 0]]) +
          test.Cubic('Phiy', eval=test.Samples(phiyAB, 0, l, 5)) +
          test.Quartic('Uz', eval=test.Samples(uzAB, 0, l, 5)),
      BC: test.Linear('Vz', -q * l / 2, -3 / 2 * q * l) +
          test.Quadratic('My', eval=[[0, 0], [l, -q * l * l]]) +
          test.Cubic('Phiy', eval=test.Samples(phiyBC, 0, l, 5)) +
          test.Quartic('Uz', eval=test.Samples(uzBC, 0, l, 5)),
    },
  },
};

[
  beams_with_shear_hinge(q=1.2e3, l=2, E=210000e6, Iyy=10e-6, hinge='left'),
  beams_with_shear_hinge(q=1.2e3, l=2, E=210000e6, Iyy=10e-6, hinge='right'),
  beams_with_shear_hinge(q=-2e3, l=5, E=100000e6, Iyy=20e-6, hinge='left'),
  beams_with_shear_hinge(q=-2e3, l=5, E=100000e6, Iyy=20e-6, hinge='right'),
  beams_with_rot_hinge(q=1.2e3, l=2, E=210000e6, Iyy=10e-6, hinge='left'),
  beams_with_rot_hinge(q=1.2e3, l=2, E=210000e6, Iyy=10e-6, hinge='right'),
  beams_with_rot_hinge(q=-2e3, l=5, E=100000e6, Iyy=20e-6, hinge='left'),
  beams_with_rot_hinge(q=-2e3, l=5, E=100000e6, Iyy=20e-6, hinge='right'),
]
