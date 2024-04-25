# Development Guide

If you read this, you might be interested in contributing; this is great! There are many ways to
contribute to this project. You can file a new issue when you run into a bug or something that is
hard to use. You can add a new boundary value problem as an automated test. Or you can fix a bug or
implement a feature that is missing. Depending on the contribution, it can make sense to open an
issue beforehand to discuss the idea. In any case, contributions are very welcome!

## Required and optional tools

Before starting, make ensure these essential tools are installed:

- [Go](https://go.dev) to compile the solver and run tests.
- [Pre-commit](https://pre-commit.com) to run linter, formatter, and other checks (see details
  below).
- [Git](https://git-scm.com), obviously.

Depending on the workflow, these optional tools might make sense, too:

- [Python](https://www.python.org) and [Sympy](https://www.sympy.org) for interpreting utilities in
  `script/`. Python version 3.12 is expected to work well for now, older versions might work, too.
- [Jsonnet](https://jsonnet.org) and [jq](https://jqlang.github.io/jq) for working with integration
  test specifications. While the test binary links against a Jsonnet package to evaluate these files
  on the fly, both tools can be helpful for inspection and debugging `.jsonnet` files directly.
- [wazero](https://wazero.io) for running Webassembly binaries. Any binary installation should work
  (package manager, `go install`, etc.). [wasmtime](https://wasmtime.dev) should generally work as
  well, it doesn't seem to have an option to mount files from the filesystem into the runtime yet -
  which is required to run the integration test suite.
- If you like `Makefile`s then `make`, e.g. [GNU Make](https://www.gnu.org/software/make) or BSD
  make, at this point both should work as the Makefile is so minimal that it is almost useless.

### Install pre-commit hooks

It's best to run the linting and formatting steps locally, to avoid CI surprises and to have a
clean, `bisect`-able history. First, install `pre-commit` itself, see the [pre-commit installation
docs](https://pre-commit.com/#install). When `pre-commit` is in your `$PATH`, make sure the
configured hooks are run automatically, e.g., using these commands from the repo root:
```shell
echo -e '#!/usr/bin/env sh\n\n./scripts/clean_staged_pre_commit.sh $@' > .git/hooks/pre-commit
chmod +x !$
```
This causes `scripts/clean_staged_pre_commit.sh` to run locally before every commit. This custom
script is different from what the upstream `pre-commit` docs suggest. The reason for that is the
conceptual difference between the file-centric `pre-commit` tool and the package-centric Go
ecosystem. Go tooling will process every file in a package, whether `git`-tracked or not, and this
can lead to false negatives or false positives when running `pre-commit` as is. The custom script
ensures that only tracked files with their staged changes are subject to the pre-commit hooks, which
mirrors exactly what CI would check.

## Building and running tests

There is nothing special about how the repository is organised. To compile the solver and related
tooling,
```shell
go build ./...
```
and to run all tests:
```shell
go test ./...
```
There's a `Makefile`, too, which shows how to invoke some variations of the above or related
commands, e.g., to report test coverage or to build and test Wasm binaries.

## Pull requests

For any patches other than trivialities, please provide a meaningful PR description that explains
what is proposed and why. It can happen that a PR is not merged after it has been reviewed. Avoiding
feature creep or maintaining the architectural consistency of the project are valid reasons to not
merge a contribution. That being said, this should be an unlikely scenario, and should certainly not
keep you from opening PRs. If unsure, open an issue first to spawn a discussion around the change or
addition you have in mind.

## Decision records

There are several [ADR](https://adr.github.io)s in [decisions/](decisions). These might might help
understanding aspects of the implementation, and it is recommended to read these to get a feeling
for how the repository works.
