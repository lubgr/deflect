local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

local default = bvp.Defaults();

local frame_inclined(l, w, A, Iyy) = {
  name: 'frame_inclined_support',
  description: 'Horizontal frame with inclined support',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Generic('default', A=A, Iyy=Iyy),

  elements: {
    AB: default.Frame2d(),
  },

  dirichlet: {
    A: bvp.Ux() + bvp.Phiy() + bvp.Uz(w),
  },

  local alpha = 45 * bvp.pi / 180.0,

  links: {
    B: bvp.InclinedSupportUxUz(alpha),
  },

  expected: {
    primary: {
      local x0 = A * l * l / (2 * Iyy),
      local x1 = w / ((2 * x0 + 3) * l),
      B: test.Ux(x1 * 3 * l) + test.Uz(x1 * 3 * l) + test.Phiy(x1 * 3 * x0),
    },
  },
};

[
  frame_inclined(l=2.0, w=0.01, A=0.01, Iyy=10e-6),
  frame_inclined(l=3.0, w=0.005, A=0.005, Iyy=25e-6),
]
