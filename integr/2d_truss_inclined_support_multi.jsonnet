local bvp = import 'bvp.libsonnet';
local test = import 'test.libsonnet';

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

  material: bvp.LinElast('default', E=30000e6, nu=0.3, rho=1),
  crosssection: bvp.Rectangle('default', b=0.1, h=0.1),

  elements: {
    AN1: bvp.Truss2d(nodes=['A', 'N1']),
    N1B: bvp.Truss2d(nodes=['N1', 'B']),
    AC: bvp.Truss2d(nodes=['A', 'C']),
    CN2: bvp.Truss2d(nodes=['C', 'N2']),
    N2N3: bvp.Truss2d(nodes=['N2', 'N3']),
    N1N2: bvp.Truss2d(nodes=['N1', 'N2']),
    AN2: bvp.Truss2d(nodes=['A', 'N2']),
    CN1: bvp.Truss2d(nodes=['C', 'N1']),
    N1N3: bvp.Truss2d(nodes=['N1', 'N3']),
    N2B: bvp.Truss2d(nodes=['N2', 'B']),
  },

  neumann: {
    N1: bvp.Fx(2 * F) + bvp.Fz(F),
    N2: bvp.Fz(-F),
    N3: bvp.Fx(F) + bvp.Fz(-2 * F),
  },

  links: {
    A: bvp.InclinedSupportUxUz(alpha * bvp.pi / 180.0),
    B: bvp.InclinedSupportUxUz(beta * bvp.pi / 180.0),
    C: bvp.InclinedSupportUxUz(gamma * bvp.pi / 180.0),
  },
};

// The reaction forces for variables angles alpha, beta, and gamma are the solution to this system
// of equations, where sa = sin(alpha), ca = cos(alpha), sb = sin(beta), ...:
//   A·sa + B·sb   + C·sg = 3·F
//   A·ca + B·cb   + C·cg = 2·F
//          B·2·cb + C·sg = 5·F
// The following Python/numpy snippet solves for A, B, and C, and decomposes them into x/z
// components:
//   A = np.array([[sa,  sb, sg], [ca,  cb, cg], [0, 2*cb, sg]])
//   b = np.array([3*F, 2*F, 5*F])
//   x = np.linalg.solve(A, b)
//   A = (-x[0]*sa, x[0]*ca)
//   B = (-x[1]*sb, x[1]*cb)
//   C = (-x[2]*sg, x[2]*cg)
//   print(f"  A: test.Fx({A[0]}) + test.Fz({A[1]}),")
//   print(f"  B: test.Fx({B[0]}) + test.Fz({B[1]}),")
//   print(f"  C: test.Fx({C[0]}) + test.Fz({C[1]}),")

local bvp0 = structure(alpha=45, beta=100, gamma=90, a=1) + {
  expected: {
    reaction: {
      A: test.Fx(-15387.07185026498) + test.Fz(15387.071850264982),
      B: test.Fx(26161.215550794954) + test.Fz(4612.928149735015),
      C: test.Fx(-40774.143700529974) + test.Fz(2.4966962285414125e-12),
    },
  },
};

local bvp1 = structure(alpha=0, beta=0, gamma=90, a=0.9876) + {
  expected: {
    reaction: {
      A: test.Fx(0) + test.Fz(10e3),
      B: test.Fx(0) + test.Fz(10e3),
      C: test.Fx(-30e3) + test.Fz(0),
    },
  },
};

local bvp2 = structure(alpha=90, beta=0, gamma=90, a=10) + {
  expected: {
    reaction: {
      A: test.Fx(-20e3) + test.Fz(0),
      B: test.Fx(0) + test.Fz(20e3),
      C: test.Fx(-10e3) + test.Fz(0),
    },
  },
};

local bvp3 = structure(alpha=15, beta=25, gamma=35, a=1) + {
  expected: {
    reaction: {
      A: test.Fx(10786.549852769092) + test.Fz(-40255.952088908845),
      B: test.Fx(-2801.2804422787344) + test.Fz(6007.365294754823),
      C: test.Fx(-37985.26941049036) + test.Fz(54248.58679415402),
    },
  },
};

local bvp4 = structure(alpha=200, beta=250, gamma=300, a=1) + {
  expected: {
    reaction: {
      A: test.Fx(788149.3380187545) + test.Fz(-2165422.509364637),
      B: test.Fx(-2823460.4899709313) + test.Fz(1027655.5759760885),
      C: test.Fx(2005311.1519521773) + test.Fz(1157766.9333885484),
    },
  },
};

[bvp0, bvp1, bvp2, bvp3, bvp4]
