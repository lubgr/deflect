{
  pi:: std.acos(0) * 2,

  Ux(value=0):: [{ Ux: value }],
  Uy(value=0):: [{ Uy: value }],
  Uz(value=0):: [{ Uz: value }],
  Phix(value=0):: [{ Phix: value }],
  Phiy(value=0):: [{ Phiy: value }],
  Phiz(value=0):: [{ Phiz: value }],

  InclinedSupportUxUz(angle):: [{ from: 'Ux', to: 'Uz', angle: angle }],
  InclinedSupportUxUy(angle):: [{ from: 'Ux', to: 'Uy', angle: angle }],
  InclinedSupportUyUz(angle):: [{ from: 'Uy', to: 'Uz', angle: angle }],

  local single(what, value, x) =
    if x == null then
      {
        [what]: value,
      }
    else
      {
        kind: what,
        degree: 'dirac',
        position: [x],
        values: [value],
      },

  Fx(value, x=null):: [single('Fx', value, x)],
  Fy(value, x=null):: [single('Fy', value, x)],
  Fz(value, x=null):: [single('Fz', value, x)],
  Mx(value, x=null):: [single('Mx', value, x)],
  My(value, x=null):: [single('My', value, x)],
  Mz(value, x=null):: [single('Mz', value, x)],

  local constant(what, value) =
    assert std.isNumber(value);
    {
      kind: what,
      degree: 'constant',
      values: [value],
    },

  local linear(what, values, x) =
    assert std.isArray(values);
    {
      kind: what,
      degree: 'linear',
      values: values,
      [if x != null then 'position']: x,
    },

  local dispatch(what, values, x) =
    if std.isNumber(values) then
      constant(what, values)
    else
      linear(what, values, x),

  qx(values, x=null):: [dispatch('qx', values, x)],
  qy(values, x=null):: [dispatch('qy', values, x)],
  qz(values, x=null):: [dispatch('qz', values, x)],
  mx(values, x=null):: [dispatch('mx', values, x)],
  my(values, x=null):: [dispatch('my', values, x)],
  mz(values, x=null):: [dispatch('mz', values, x)],

  LinElast(id, E, nu, rho)::
    {
      [id]: {
        kind: 'linelast',
        parameter: {
          E: E,
          nu: nu,
          rho: rho,
        },
      },
    },

  Rectangle(id, b, h)::
    {
      [id]: {
        kind: 'rectangle',
        parameter: {
          b: b,
          h: h,
        },
      },
    },

  Generic(id, A, Iyy)::
    {
      [id]: {
        kind: 'constants',
        parameter: {
          A: A,
          Iyy: Iyy,
        },
      },
    },

  Defaults(material='default', cs='default')::
    {
      local element(nodes, hinges, kind) = {
        kind: kind,
        material: material,
        cs: cs,
        [if std.length(nodes) > 0 then 'nodes']: nodes,
        [if std.length(hinges) > 0 then 'hinges']: hinges,
      },

      Truss2d(nodes=[], hinges={})::
        element(nodes, hinges, 'truss2d'),
      Truss3d(nodes=[], hinges={})::
        element(nodes, hinges, 'truss3d'),
      Frame2d(nodes=[], hinges={})::
        element(nodes, hinges, 'frame2d'),
    },

  local allowedPolynomials = ['Ux', 'Uz', 'Uy', 'Phiy', 'Phiz', 'Phix', 'Nx', 'Vz', 'Vy', 'My', 'Mz', 'Mx'],

  Constant(kind, value, range=null)::
    assert std.member(allowedPolynomials, kind) : "Unknown function '%s'" % kind;
    [
      {
        kind: kind,
        degree: 0,
        boundary: [value],
        [if range != null then 'range']: range,
      },
    ],

  Linear(kind, first, second, range=null)::
    assert std.member(allowedPolynomials, kind) : "Unknown function '%s'" % kind;
    [
      {
        kind: kind,
        degree: 1,
        boundary: [first, second],
        [if range != null then 'range']: range,
      },
    ],

  local higherOrder(kind, degree, eval, range) = {
    assert std.member(allowedPolynomials, kind) : "Unknown function '%s'" % kind,
    kind: kind,
    degree: degree,
    [if range != null then 'range']: range,
    [if eval != null then 'eval']: eval,
  },

  Quadratic(kind, eval=null, range=null):: [higherOrder(kind, 2, eval, range)],
  Cubic(kind, eval=null, range=null):: [higherOrder(kind, 3, eval, range)],
  Quartic(kind, eval=null, range=null):: [higherOrder(kind, 4, eval, range)],
  Quintic(kind, eval=null, range=null):: [higherOrder(kind, 5, eval, range)],

  Samples(f, x0, xE, n):: [
    [x, f(x)]
    for x in [x0 + i * (xE - x0) / n for i in std.range(0, n)]
  ],
}
