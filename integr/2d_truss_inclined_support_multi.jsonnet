local lib = import 'common.libsonnet';

local default = lib.Defaults();

local F = 10e3;

local structure(alpha, beta, gamma, a) = {
  name: 'inclined_support_multi_' + std.join('_', ['%.0f' % f for f in [alpha, beta, gamma, a]]),
  description: 'Truss structure with three inclined supports',

  nodes: {
    A: [0, 0, 0],
    B: [2 * a, 0, 0],
    C: [0, 0, a],
    N1: [a, 0, 0],
    N2: [a, 0, a],
    N3: [2 * a, 0, a],
  },

  material: lib.LinElast('default', E=30000 * 1e6, nu=0.3, rho=1),
  crosssection: lib.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AN1: default.Truss2d(nodes=['A', 'N1']),
    N1B: default.Truss2d(nodes=['N1', 'B']),
    AC: default.Truss2d(nodes=['A', 'C']),
    CN2: default.Truss2d(nodes=['C', 'N2']),
    N2N3: default.Truss2d(nodes=['N2', 'N3']),
    N1N2: default.Truss2d(nodes=['N1', 'N2']),
    AN2: default.Truss2d(nodes=['A', 'N2']),
    CN1: default.Truss2d(nodes=['C', 'N1']),
    N1N3: default.Truss2d(nodes=['N1', 'N3']),
    N2B: default.Truss2d(nodes=['N2', 'B']),
  },

  neumann: {
    N1: lib.Fx(2 * F) + lib.Fz(F),
    N2: lib.Fz(-F),
    N3: lib.Fx(F) + lib.Fz(-2 * F),
  },

  links: {
    A: lib.LinkUxUzAngular(alpha * lib.pi / 180.0),
    B: lib.LinkUxUzAngular(beta * lib.pi / 180.0),
    C: lib.LinkUxUzAngular(gamma * lib.pi / 180.0),
  },
};

// The reaction forces for variables angles alpha, beta, and gamma are the solution to this system
// of equations, where sa = sin(alpha), ca = cos(alpha), sb = sin(beta), ...:
//   A*sa + B*sb   + C*sg = 3*F
//   A*ca + B*cb   + C*cg = 2*F
//          B*2*cb + C*sg = 5*F
// The following Python/numpy snippet solves for A, B, and C, and decomposes them into x/z
// components:
//   A = np.array([[sa,  sb, sg], [ca,  cb, cg], [0, 2*cb, sg]])
//   b = np.array([3*F, 2*F, 5*F])
//   x = np.linalg.solve(A, b)
//   A = (-x[0]*sa, x[0]*ca)
//   B = (-x[1]*sb, x[1]*cb)
//   C = (-x[2]*sg, x[2]*cg)
//   print(f"  A: lib.Fx({A[0]}) + lib.Fz({A[1]}),")
//   print(f"  B: lib.Fx({B[0]}) + lib.Fz({B[1]}),")
//   print(f"  C: lib.Fx({C[0]}) + lib.Fz({C[1]}),")

local bvp0 = structure(alpha=45, beta=100, gamma=90, a=1) + {
  expected: {
    reaction: {
      A: lib.Fx(-15387.07185026498) + lib.Fz(15387.071850264982),
      B: lib.Fx(26161.215550794954) + lib.Fz(4612.928149735015),
      C: lib.Fx(-40774.143700529974) + lib.Fz(2.4966962285414125e-12),
    },
  },
};

local bvp1 = structure(alpha=0, beta=0, gamma=90, a=0.9876) + {
  expected: {
    reaction: {
      A: lib.Fx(0) + lib.Fz(10e3),
      B: lib.Fx(0) + lib.Fz(10e3),
      C: lib.Fx(-30e3) + lib.Fz(0),
    },
  },
};

local bvp2 = structure(alpha=90, beta=0, gamma=90, a=10) + {
  expected: {
    reaction: {
      A: lib.Fx(-20e3) + lib.Fz(0),
      B: lib.Fx(0) + lib.Fz(20e3),
      C: lib.Fx(-10e3) + lib.Fz(0),
    },
  },
};

local bvp3 = structure(alpha=15, beta=25, gamma=35, a=1) + {
  expected: {
    reaction: {
      A: lib.Fx(10786.549852769092) + lib.Fz(-40255.952088908845),
      B: lib.Fx(-2801.2804422787344) + lib.Fz(6007.365294754823),
      C: lib.Fx(-37985.26941049036) + lib.Fz(54248.58679415402),
    },
  },
};

local bvp4 = structure(alpha=200, beta=250, gamma=300, a=1) + {
  expected: {
    reaction: {
      A: lib.Fx(788149.3380187545) + lib.Fz(-2165422.509364637),
      B: lib.Fx(-2823460.4899709313) + lib.Fz(1027655.5759760885),
      C: lib.Fx(2005311.1519521773) + lib.Fz(1157766.9333885484),
    },
  },
};

[bvp0, bvp1, bvp2, bvp3, bvp4]
