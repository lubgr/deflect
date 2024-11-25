#!/usr/bin/env python3

import argparse
import sys
import textwrap

from sympy import Matrix, Symbol, ccode, eye, pretty, symbols, zeros
from sympy.codegen.rewriting import create_expand_pow_optimization


def show(description, t, k, ksym, r, d):
    indent = " " * 4
    print(f"{description}, tangent in global coordinates:\n")
    if k != None:
        print(textwrap.indent(pretty(t.T * k * t), indent))
        print()
    print(
        textwrap.indent(
            pretty(t.T * ksym * t, wrap_line=False, use_unicode=False), indent
        )
    )
    print(f"\n{description}, rhs in global coordinates:\n")
    print(textwrap.indent(pretty(t.T * r), indent))
    print()
    print(f"\n{description}, d in local coordinates:\n")
    print(textwrap.indent(pretty(t * d), indent))
    print()


def code(description, t, ksym, r, d, indices):
    indent = " " * 4
    expandpow = create_expand_pow_optimization(6)
    kglobal = t.T * ksym * t
    rglobal = t.T * r
    size = kglobal.rows

    print(f"{description}, tangent in global coordinates:\n")

    for i in range(0, size):
        isym = indices[i]
        for j in range(i, size):
            jsym = indices[j]
            code = ccode(expandpow(kglobal[i, j]))

            print(textwrap.indent(f"({isym}, {jsym}) = {code}", indent))

    print(f"\n{description}, rhs in global coordinates:\n")

    for i in range(0, size):
        isym = indices[i]
        code = ccode(expandpow(rglobal[i]))
        print(textwrap.indent(f"({isym}) = {code}", indent))


def truss3d():
    EA, l, cx, cy, cz = symbols("EA l cx cy cz")

    k = zeros(2, 2)
    k[0, 0] = k[1, 1] = EA / l
    k[0, 1] = k[1, 0] = -EA / l

    t = zeros(2, 6)
    t[0, 0] = t[1, 3] = cx
    t[0, 1] = t[1, 4] = cy
    t[0, 2] = t[1, 5] = cz

    ksym = Matrix(2, 2, lambda i, j: Symbol(f"k{i}{j}"))
    r = Matrix(2, 1, lambda i, _: Symbol(f"r{i}"))
    d = Matrix(6, 1, lambda i, _: Symbol(f"d{i}"))

    show("3d truss", t, k, ksym, r, d)


def truss2d():
    EA, l, c, s = symbols("EA l c s")
    size = 4

    k = zeros(size, size)
    k[0, 0] = k[2, 2] = EA / l
    k[0, 2] = k[2, 0] = -EA / l

    ksym = zeros(size, size)
    ksym[0, 0] = Symbol("k00")
    ksym[0, 2] = Symbol("k01")
    ksym[2, 0] = Symbol("k01")
    ksym[2, 2] = Symbol("k11")

    t = eye(size, size)
    t[0, 0] = t[2, 2] = c
    t[0, 1] = t[2, 3] = s
    t[1, 0] = t[3, 2] = s
    t[1, 1] = t[3, 3] = -c

    r = zeros(size, 1)
    r[0] = Symbol("r1")
    r[1] = 0
    r[2] = Symbol("r2")
    r[3] = 0

    d = Matrix(size, 1, lambda i, _: Symbol(f"d{i}"))

    show("2d truss", t, k, ksym, r, d)


def beam3d():
    kh = lambda i, j: Symbol(f"kh({i},{j})")
    kv = lambda i, j: Symbol(f"kv({i},{j})")
    rh = lambda i: Symbol(f"rh({i})")
    rv = lambda i: Symbol(f"rv({i})")

    size = 12
    ksym = zeros(size, size)

    ksym[1, 1] = kh(0, 0)
    ksym[1, 5] = kh(0, 1)
    ksym[1, 7] = kh(0, 2)
    ksym[1, 11] = kh(0, 3)
    ksym[5, 5] = kh(1, 1)
    ksym[5, 7] = kh(1, 2)
    ksym[5, 11] = kh(1, 3)
    ksym[7, 7] = kh(2, 2)
    ksym[7, 11] = kh(2, 3)
    ksym[11, 11] = kh(3, 3)

    ksym[2, 2] = kv(0, 0)
    ksym[2, 4] = kv(0, 1)
    ksym[2, 8] = kv(0, 2)
    ksym[2, 10] = kv(0, 3)
    ksym[4, 4] = kv(1, 1)
    ksym[4, 8] = kv(1, 2)
    ksym[4, 10] = kv(1, 3)
    ksym[8, 8] = kv(2, 2)
    ksym[8, 10] = kv(2, 3)
    ksym[10, 10] = kv(3, 3)

    # Make matrix symmetric
    for i in range(0, size):
        for j in range(0, i):
            ksym[i, j] = ksym[j, i]

    ts = Matrix(3, 3, lambda i, j: Symbol(f"t{i}{j}"))
    t = zeros(size, size)
    t[0:3, 0:3] = ts
    t[3:6, 3:6] = ts
    t[6:9, 6:9] = ts
    t[9:12, 9:12] = ts

    r = zeros(size, 1)
    r[1] = rh(0)
    r[2] = rv(0)
    r[4] = rv(1)
    r[5] = rh(1)
    r[7] = rh(2)
    r[8] = rv(2)
    r[10] = rv(3)
    r[11] = rh(3)

    d = Matrix(size, 1, lambda i, _: Symbol(f"d{i}"))
    indices = zeros(size, 1)

    indices[0] = Symbol("ux0")
    indices[1] = Symbol("uy0")
    indices[2] = Symbol("uz0")
    indices[3] = Symbol("phix0")
    indices[4] = Symbol("phiy0")
    indices[5] = Symbol("phiz0")
    indices[6] = Symbol("ux1")
    indices[7] = Symbol("uy1")
    indices[8] = Symbol("uz1")
    indices[9] = Symbol("phix1")
    indices[10] = Symbol("phiy1")
    indices[11] = Symbol("phiz1")

    show("3d beam", t, None, ksym, r, d)
    code("3d beam", t, ksym, r, d, indices)


def beam2d():
    size = 6
    EI, l, c, s = symbols("EI l c s")
    l2, l3 = l * l, l * l * l

    k = zeros(size, size)
    k[1, 1] = 12 * EI / l3
    k[1, 2] = k[2, 1] = -6 * EI / l2
    k[1, 4] = k[4, 1] = -12 * EI / l3
    k[1, 5] = k[5, 1] = -6 * EI / l2
    k[2, 2] = 4 * EI / l
    k[2, 4] = k[4, 2] = 6 * EI / l2
    k[2, 5] = k[5, 2] = 2 * EI / l
    k[4, 4] = 12 * EI / l3
    k[4, 5] = k[5, 4] = 6 * EI / l2
    k[5, 5] = 4 * EI / l

    ksym = zeros(size, size)
    ksym[1, 1] = Symbol("k00")
    ksym[1, 2] = ksym[2, 1] = Symbol("k01")
    ksym[1, 4] = ksym[4, 1] = Symbol("k02")
    ksym[1, 5] = ksym[5, 1] = Symbol("k03")
    ksym[2, 2] = Symbol("k11")
    ksym[2, 4] = ksym[4, 2] = Symbol("k12")
    ksym[2, 5] = ksym[5, 2] = Symbol("k13")
    ksym[4, 4] = Symbol("k22")
    ksym[4, 5] = ksym[5, 4] = Symbol("k23")
    ksym[5, 5] = Symbol("k33")

    r = zeros(size, 1)
    r[0] = 0
    r[1] = Symbol("r1")
    r[2] = Symbol("r2")
    r[3] = 0
    r[4] = Symbol("r4")
    r[5] = Symbol("r5")

    d = Matrix(size, 1, lambda i, _: Symbol(f"d{i}"))

    t = eye(size, size)
    t[0, 0] = t[3, 3] = c
    t[0, 1] = t[3, 4] = s
    t[1, 0] = t[4, 3] = s
    t[1, 1] = t[4, 4] = -c
    t[2, 2] = t[5, 5] = -1

    show("2d beam", t, k, ksym, r, d)


parser = argparse.ArgumentParser(
    description="""
    Prints element matrices transformed from the local element CO system to the global
    one. This can help implementing fast element implementations, since using the
    resulting matrix expressions saves matrix-matrix and matrix-vector products with
    complete transformation matrices.
    """
)
parser.add_argument(
    "--element",
    type=str,
    choices=["beam2d", "beam3d", "truss2d", "truss3d"],
    required=True,
    help="Show results for this element type",
)
args = parser.parse_args()

fct = getattr(sys.modules[__name__], args.element)
fct()
