# Implement the Solver in Go

## Context and Problem Statement

Choose a language for implementing element interpolation, matrix assembly, boundary condition
handling and an equation solver.

## Considered Options

* Go
* C++
* Zig
* Ocaml
* Rust

## Decision Outcome

This is a hobby project, so there is no actual pressure to justify any language choice.

However, how fun the project is long term does depend on sound technological choices. It should be
easy to pick up after some time or to contribute a patch without being a primary maintainer. We also
want a portable solution with support for WebAssembly, since a deployment with web technologies is
an option for the future. Performance is important (when is it not), but not to a degree where any
complexity is acceptable to squeeze out whatever possible.

Golang seems to fit nicely here. It fulfils the criteria above and is known for its gentle learning
curve. While Zig has many similar characteristics, it is not as mature at the time of writing. Ocaml
is very mature, but less common, which could keep contributors from providing patches. Given the
performance requirements, the complexity of Rust or C++ are not justified.

Again, any decision would be okay, and this file exists merely to document some of the initial ideas
around the choice.
