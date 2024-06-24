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

  Rectangle(id, b, h, roll=0)::
    {
      [id]: {
        kind: 'rectangle',
        parameter: {
          b: b,
          h: h,
          [if roll != 0 then 'roll']: roll,
        },
      },
    },

  Generic(id, A, Iyy, Izz, roll=0)::
    {
      [id]: {
        kind: 'constants',
        parameter: {
          A: A,
          Iyy: Iyy,
          Izz: Izz,
          [if roll != 0 then 'roll']: roll,
        },
      },
    },

  local element(nodes, hinges, kind, material, cs) = {
    kind: kind,
    [if material != null then 'material']: material,
    [if cs != null then 'cs']: cs,
    [if std.length(nodes) > 0 then 'nodes']: nodes,
    [if std.length(hinges) > 0 then 'hinges']: hinges,
  },

  Truss2d(nodes=[], hinges={}, material=null, cs=null)::
    element(nodes, hinges, 'truss2d', material, cs),
  Truss3d(nodes=[], hinges={}, material=null, cs=null)::
    element(nodes, hinges, 'truss3d', material, cs),
  Frame2d(nodes=[], hinges={}, material=null, cs=null)::
    element(nodes, hinges, 'frame2d', material, cs),
}
