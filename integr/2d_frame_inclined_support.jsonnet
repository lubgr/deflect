local lib = import 'common.libsonnet';

local default = lib.Defaults();

local frame_inclined(l, w, A, Iyy) = {
  name: 'frame_inclined_support',
  description: 'Horizontal frame with inclined support',

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Generic('default', A=A, Iyy=Iyy),

  elements: {
    AB: default.Frame2d(),
  },

  dirichlet: {
    A: lib.Ux() + lib.Phiy() + lib.Uz(w),
  },

  local alpha = 45 * lib.pi / 180.0,

  links: {
    B: lib.InclinedSupportUxUz(alpha),
  },

  expected: {
    primary: {
      local x0 = A * l * l / (2 * Iyy),
      local x1 = w / ((2 * x0 + 3) * l),
      B: lib.Ux(x1 * 3 * l) + lib.Uz(x1 * 3 * l) + lib.Phiy(x1 * 3 * x0),
    },
  },
};

[
  frame_inclined(l=2.0, w=0.01, A=0.01, Iyy=10e-6),
  frame_inclined(l=3.0, w=0.005, A=0.005, Iyy=25e-6),
]
