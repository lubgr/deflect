# How to Specify Data for Integration/End-to-End Tests?

## Context and Problem Statement

Simple boundary value problems are natural test cases for this project. Ideally, we can easily add
new tests and maintain existing ones, without being held up by verbose syntax, a restrictive data
format, or other factors. One or more test cases could be specified in a single file, listing input
data (nodes, material and cross sections, elements, boundary conditions) and expected results (nodal
results, interpolations).

## Considered Options

- [Dhall](https://dhall-lang.org)
- Standard table-driven Go tests with native helper APIs
- Home grown plain text format
- [CUE](https://cuelang.org)
- [Jsonnet](https://jsonnet.org)
- Bindings for a scripting language, e.g. [Gopher-Lua](https://github.com/yuin/gopher-lua)

## Decision Outcome

We are going with Jsonnet. It is close to JSON, which could be a data exchange format for a UI. This
will allow us to easily use the test data for ordinary integration tests as well as for tests closer
to a deployment, e.g. when a service is to expose a REST API. Jsonnet has straightforward Go
integration, because the most recent official implementation is written in Go. The flexibility of
the language offers enough room for DRY. Arguably, it's too flexible for a configuration language
(inheritance, mixins, syntactic sugar, ...), so that large configuration systems can become hard to
maintain. In this case though, we shall have short files with small amounts of shared context.

As for the others:

- `Dhall` doesn't have proper maths support, e.g. trigonometric functions are missing, and we might
  need them. Also, the Go bindings seem unmaintained. It's possible to not use the bindings at all
  and convert everything to JSON or Yaml before reading it in, but then we depend on an external
  tool to be present.
- Standard table-driven Go tests are too verbose here, and it would not be easy enough to quickly
  add new tests.
- A home grown plain text format can work, but is at risk of being either too restrictive or
  suffering from feature creep, such that over time, we duplicate what a distinct configuration
  language is intended to solve in the first place. For example, would we support variable
  definitions and arithmetic to easily specify mesh coordinates?
- `CUE` looks powerful and offers natural Go integration. However, it is so broadly applicable that
  the docs are hard to understand, and the focus seems to be missing. Not KISS enough, setting it up
  properly might not be a bit of a chore.
- Finally, binding a scripting language can be fun, and there is good malleability since the binding
  mechanics are controlled directly. However, the overhead can be substantial, and only a small
  portion of the language power would actually be in use.
