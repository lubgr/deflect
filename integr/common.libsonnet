{
  pi:: std.acos(0) * 2,

  Ux(value=0):: [{ Ux: value }],
  Uy(value=0):: [{ Uy: value }],
  Uz(value=0):: [{ Uz: value }],
  Phix(value=0):: [{ Phix: value }],
  Phiy(value=0):: [{ Phiy: value }],
  Phiz(value=0):: [{ Phiz: value }],

  LinkUxUzAngular(angle):: [{ from: 'Ux', to: 'Uz', angle: angle }],
  LinkUxUyAngular(angle):: [{ from: 'Ux', to: 'Uy', angle: angle }],
  LinkUyUzAngular(angle):: [{ from: 'Uy', to: 'Uz', angle: angle }],

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

  Defaults(material='default', cs='default')::
    {
      local common = {
        material: material,
        cs: cs,
      },

      Truss2d(nodes=[], disjoint=[])::
        common {
          kind: 'truss2d',
          [if std.length(nodes) > 0 then 'nodes']: nodes,
          [if std.length(disjoint) > 0 then 'disjoint']: disjoint,
        },
    },
}
