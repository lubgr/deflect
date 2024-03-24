#!/usr/bin/env python3

import argparse
import copy
import textwrap

from sympy import Matrix, Symbol, ccode, pretty, simplify, zeros
from sympy.codegen.rewriting import create_expand_pow_optimization


def create_matrices(size):
    k, r, d = zeros(size, size), zeros(size, 1), zeros(size, 1)

    for i in range(0, size):
        r[i] = Symbol(f"r{i}")
        d[i] = Symbol(f"d{i}")
        for j in range(i, size):
            # It would be more concise to construct k with a MatrixSymbol, but that
            # leads to undesired output when generating code, where we want plain dumb
            # symbols like 'k23', not sympy-magical matrix comprehension logic.
            k[i, j] = k[j, i] = Symbol(f"k{i}{j}")

    return k, d, r


def layout(size, hinges):
    assert sorted(hinges) == hinges
    result = copy.deepcopy(hinges)

    for i in range(0, size):
        if i not in hinges:
            result.append(i)

    return result


def reorder(k, d, r, hinges):
    size = k.rows
    tmpk1, tmpk2 = zeros(size, size), zeros(size, size)
    tmpd, tmpr = zeros(size, 1), zeros(size, 1)
    reorder = layout(size, hinges)

    for i, source in enumerate(reorder):
        tmpk1[i, :] = k[source, :]
        tmpd[i, :] = d[source, :]
        tmpr[i, :] = r[source, :]
    for i, source in enumerate(reorder):
        tmpk2[:, i] = tmpk1[:, source]

    return tmpk2, tmpd, tmpr


def show(k, r, d, hinges):
    print(f"Size: {k.rows}, hinges: {hinges}\n")

    def indent(txt):
        print(textwrap.indent(txt, 4 * " "))

    indent(pretty(k))
    print()

    expandpow = create_expand_pow_optimization(6)
    for i in range(0, k.rows):
        for j in range(i, k.cols):
            code = ccode(expandpow(k[i, j]))
            indent(f"k.SetSym({i}, {j}, {code})")

    for name, vec in [("r", r), ("d", d)]:
        print()
        indent(pretty(vec))
        print()

        for i in range(0, vec.rows):
            code = ccode(expandpow(vec[i]))
            indent(f"{name}.SetVec({i}, {code})")


def condensate(size, hinges):
    k, d, r = create_matrices(size)
    k, d, r = reorder(k, d, r, hinges)

    n = len(hinges)
    k11, k22 = Matrix(k[:n, :n]), Matrix(k[n:, n:])
    k12, k21 = Matrix(k[:n, n:]), Matrix(k[n:, :n])
    r1, r2 = Matrix(r[:n]), Matrix(r[n:])
    d1, d2 = Matrix(d[:n]), Matrix(d[n:])

    k11inv = k11.inv()
    kred = k22 - k21 * k11inv * k12
    rred = r2 - k21 * k11inv * r1
    d1 = k11inv * (r1 - k12 * d2)

    kred, rred, d1 = simplify(kred), simplify(rred), simplify(d1)
    _, d, _ = create_matrices(size)

    for r, i in enumerate(hinges):
        cols = kred.cols
        kred = kred.row_insert(i, zeros(1, cols))
        rows = kred.rows
        kred = kred.col_insert(i, zeros(rows, 1))
        rred = rred.row_insert(i, zeros(1, 1))
        d[i] = d1[r]

    show(kred, rred, d, hinges)


parser = argparse.ArgumentParser(
    description="""
    Computes element matrices for elements with hinges, using static condensation.
    Prints the matrices along with example code. The resulting matrices are expressed in
    terms of the original element matrices for the case without any hinges, so that an
    element implementation can keep static condensation separate from the rest (material
    law, interpolation degree). Static condensation here refers to a partitioning of the
    system k·d = r into two subsystems, and eliminating the degrees of freedom that are
    hinged. For more info, see e.g. Chapter 8.1 in 'Concepts and Applications of Finite
    Element Analysis', 3rd ed., by Cook, Malkus & Plesha.
"""
)
parser.add_argument(
    "--size",
    type=int,
    required=True,
    help="""The size of the system of equations, e.g., 2 for a two-dimensional truss, 4
    for a two-dimensional beam.""",
)
parser.add_argument(
    "--hinges",
    nargs="+",
    type=int,
    required=True,
    help="""Matrix indices of the hinges that are statically condensed out. Must be less
    than the given size.""",
)
args = parser.parse_args()
condensate(args.size, args.hinges)
