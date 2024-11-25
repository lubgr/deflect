local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local l = 1.0;

local bar = {
  name: 'single_torsion_bar',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  elements: {
    AB: bvp.TorsionBar(),
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('beam', A=0.01, Iyy=1e-5, Izz=1e-5, Ixx=1e-5),

  dirichlet: {
    A: bvp.Phix(),
  },

  neumann: {
    B: bvp.Mx(0.1),
  },

  expected: {
    reaction: {
      A: test.Mx(-0.1),
    },
  },
};

[
  bar,
]
