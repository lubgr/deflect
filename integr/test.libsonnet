local bvp = import 'bvp.libsonnet';

{
  // Currently implemented identical to the bvp library, but in case this changes we keep them in
  // different libraries.
  Ux(value=0):: bvp.Ux(value),
  Uy(value=0):: bvp.Uy(value),
  Uz(value=0):: bvp.Uz(value),
  Phix(value=0):: bvp.Phix(value),
  Phiy(value=0):: bvp.Phiy(value),
  Phiz(value=0):: bvp.Phiz(value),
  // Same as the primary primitives above, implemented in terms of the bvp library.
  Fx(value, x=null):: bvp.Fx(value, x),
  Fy(value, x=null):: bvp.Fy(value, x),
  Fz(value, x=null):: bvp.Fz(value, x),
  Mx(value, x=null):: bvp.Mx(value, x),
  My(value, x=null):: bvp.My(value, x),
  Mz(value, x=null):: bvp.Mz(value, x),

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
