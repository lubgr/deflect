#!/usr/bin/env python3

import argparse
import sys
import textwrap

from sympy import Matrix, Symbol, eye, pretty, zeros


def show(description, t, k, ksym, r, d):
    indent = " " * 4
    print(f"{description}, tangent in global coordinates:\n")
    print(textwrap.indent(pretty(t.T * k * t), indent))
    print()
    print(textwrap.indent(pretty(t.T * ksym * t), indent))
    print(f"\n{description}, rhs in global coordinates:\n")
    print(textwrap.indent(pretty(t.T * r), indent))
    print()
    print(f"\n{description}, d in local coordinates:\n")
    print(textwrap.indent(pretty(t * d), indent))
    print()


def truss2d():
    size = 4
    c, s = Symbol("c"), Symbol("s")
    EA, l = Symbol("EA"), Symbol("l")

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

    d = Matrix(size, 1, lambda i, j: Symbol(f"d{i}"))

    show("2d truss", t, k, ksym, r, d)


def beam2d():
    size = 6
    c, s = Symbol("c"), Symbol("s")
    EI, l = Symbol("EI"), Symbol("l")
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

    d = Matrix(size, 1, lambda i, j: Symbol(f"d{i}"))

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
    choices=["beam2d", "truss2d"],
    required=True,
    help="Show results for this element type",
)
args = parser.parse_args()

fct = getattr(sys.modules[__name__], args.element)
fct()
