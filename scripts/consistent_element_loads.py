#!/usr/bin/env python3

import argparse
import sys
import textwrap

from sympy import (
    DiracDelta,
    Heaviside,
    Matrix,
    Symbol,
    ccode,
    diff,
    integrate,
    pretty,
    simplify,
)
from sympy.codegen.rewriting import create_expand_pow_optimization


def posSymbols(names):
    return [Symbol(name, positive=True) for name in names]


x, a, b, l = posSymbols(["x", "a", "b", "l"])
F, M, q, q0, q1 = posSymbols(["F", "M", "q", "q0", "q1"])


def show(description, matrix):
    indent = 4 * " "
    print(description + "\n")
    print(textwrap.indent(pretty(matrix), indent))

    print()

    expandpow = create_expand_pow_optimization(6)
    for i in range(0, matrix.rows):
        code = ccode(expandpow(matrix[i]))
        print(textwrap.indent(f"r{i} = {code}", indent))

    print("\n")


def distributed_loads(elmt, N):
    constant = integrate(N * q, (x, 0, l))
    show(f"Constant {elmt} load {q}:", constant)

    constant_a_b = simplify(
        integrate(N * 0, (x, 0, a))
        + integrate(N * q, (x, a, b))
        + integrate(N * 0, (x, b, l))
    )
    show(f"Constant {elmt} load {q} from {a} to {b}:", constant_a_b)

    linear = simplify(integrate(N * (q0 + (q1 - q0) * x / l), (x, 0, l)))
    show(f"Linear {elmt} load {q0}, {q1} from 0 to l:", linear)

    linear_a_b = simplify(
        integrate(N * 0, (x, 0, a))
        + integrate(N * (q0 + (q1 - q0) * x / l), (x, a, b))
        + integrate(N * 0, (x, b, l))
    )
    show(f"Linear {elmt} load {q0}, {q1} from {a} to {b}:", linear_a_b)

    triangle_l_half = simplify(
        integrate(N * q * x / (l / 2), (x, 0, l / 2))
        + integrate(N * q * (1 - (x - l / 2) / (l / 2)), (x, l / 2, l))
    )
    show(f"Linear triangular {elmt} load {q} at l/2:", triangle_l_half)

    triangle_a = simplify(
        integrate(N * q * x / a, (x, 0, a))
        + integrate(N * q * (1 - (x - a) / (l - a)), (x, a, l))
    )
    show(f"Linear triangular {elmt} load {q} at {a}:", triangle_a)


def singular_loads(elmt, N):
    f_at_a = integrate(N * DiracDelta(x - a) * F, (x, 0, l))
    f_at_a = simplify(f_at_a.subs(Heaviside(l - a), 1))
    show(f"Concentrated {elmt} force {F} at {a}:", f_at_a)


def truss2d():
    N = Matrix([(l - x) / l, x / l])
    distributed_loads("truss", N)
    singular_loads("truss", N)


def beam2d():
    N = Matrix(
        [
            1 - 3 * x**2 / l**2 + 2 * x**3 / l**3,
            -x + 2 * x**2 / l - x**3 / l**2,
            3 * x**2 / l**2 - 2 * x**3 / l**3,
            x**2 / l - x**3 / l**2,
        ]
    )

    distributed_loads("beam", N)
    singular_loads("beam", N)

    m_at_a = integrate(N * diff(DiracDelta(x - a), x) * M, (x, 0, l))
    m_at_a = simplify(m_at_a.subs(Heaviside(l - a), 1))
    m_at_a = simplify(m_at_a.subs(DiracDelta(a - l), 0))
    show(f"Concentrated beam torque {M} at {a}:", m_at_a)


parser = argparse.ArgumentParser(
    description="""
    Derives and prints example code for consistent element load vectors for common load
    cases. This shall help implementing element formulations with closed-form solutions,
    i.e., trusses, beams, or frames.
"""
)
parser.add_argument(
    "--element",
    type=str,
    choices=["beam2d", "truss2d"],
    required=True,
    help="The element type to show load vectors for",
)
args = parser.parse_args()

fct = getattr(sys.modules[__name__], args.element)
fct()
