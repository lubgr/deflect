local lib = import 'common.libsonnet';

local default = lib.Defaults();

local tetrahedon(Fz, Fy, l, b, h) = {
  // name: 'tetrahedon_%g_%g_%g_%g' % [Fz, l, b, h],
  name: 'tetrahedon_%f_%f' % [Fz, Fy],

  nodes: {
    A: [0, 0, 0],
    B: [l, 0, 0],
    C: [l, -b / 2, -h],
    D: [l, b / 2, -h],
  },

  material: lib.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AB: default.Truss3d(),
    AC: default.Truss3d(),
    AD: default.Truss3d(),
  },

  dirichlet: {
    [node]: lib.Ux() + lib.Uy() + lib.Uz()
    for node in ['B', 'C', 'D']
  },

  neumann: {
    A: lib.Fz(Fz) + lib.Fy(Fy),
  },

  expected: {
    local Rfz = {
      x: Fz * l / h,
      y: Fz * b / h / 4,
      z: Fz / 2,
    },
    local Rfy = {
      x: Fy * l / b,
      y: Fy / 2,
      z: Fy * h / b,
    },

    reaction: {
      B: lib.Fx(-Rfz.x) + lib.Fy(0) + lib.Fz(0),
      C: lib.Fx(Rfz.x / 2 + Rfy.x) + lib.Fy(-Rfz.y - Rfy.y) + lib.Fz(-Rfz.z - Rfy.z),
      D: lib.Fx(Rfz.x / 2 - Rfy.x) + lib.Fy(Rfz.y - Rfy.y) + lib.Fz(-Rfz.z + Rfy.z),
    },
    interpolation: {
      local L = std.sqrt(b * b / 4 + l * l + h * h),

      AB: lib.Constant('Nx', -Rfz.x),
      AC: lib.Constant('Nx', (Rfz.x / 2 + Rfy.x) * l / L - (-Rfz.y - Rfy.y) * b / 2 / L - (-Rfz.z - Rfy.z) * h / L),
      AD: lib.Constant('Nx', (Rfz.x / 2 - Rfy.x) * l / L + (Rfz.y - Rfy.y) * b / 2 / L - (-Rfz.z + Rfy.z) * h / L),
    },
  },
};

[
  tetrahedon(Fz=1e3, Fy=2e3, l=1, b=1, h=2),
  tetrahedon(Fz=2e3, Fy=-3e3, l=1, b=2, h=2),
  tetrahedon(Fz=-3e3, Fy=0e3, l=1, b=3, h=2),
  tetrahedon(Fz=-1e3, Fy=-2e3, l=1, b=4, h=2),
  tetrahedon(Fz=0e3, Fy=1e3, l=3, b=2, h=0.5),
  tetrahedon(Fz=1e3, Fy=2e3, l=5, b=0.132, h=0.2345),
]
