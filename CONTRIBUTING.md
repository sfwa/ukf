# Contributing to UKF

Thank you for helping. A contribution needs three things, and the maintainers
do the rest:

1. **A clear explanation** in the pull request: what changes, why, and how it
   was checked.
2. **Regression coverage** that fits the change. The unit tests run in double
   precision (`unittest`) and single precision (`unittest_single`); a change to
   filter behaviour should be tested in both. A bug fix should come with a test
   that failed before it.
3. **Documentation** for anything a user can see: the README, the
   [usage guide](docs/usage-guide.md) or the header comment, as fits.

The pull request template asks for a one-line changelog entry and what you
think the change does to compatibility. A best guess is fine. The maintainers
decide how a change is classified, which release it goes in and what version
that release gets ([RELEASING.md](RELEASING.md) says how), and write
[CHANGELOG.md](CHANGELOG.md) from the entries. You do not need to edit the
changelog, choose a version or tidy your commit history: pull requests are
usually squashed when merged.

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --target unittest unittest_single ukf_examples --parallel 2
ctest --test-dir build --output-on-failure
```

The [README](README.md#build-and-run) has the details, including how to use
another Eigen. CI builds every pull request with the compilers and Eigen
versions listed under [supported configurations](RELEASING.md#supported-configurations),
with warnings as errors, and checks the consumer fixtures in
[test/cmake](test/cmake). It is fine to open a pull request and let CI find
what your machine did not.

## Formatting and lint

Code is formatted with clang-format 18 in Google style (`.clang-format`), and
clang-tidy 18 checks it for bugs and performance problems (`.clang-tidy`). CI
runs both. Before pushing, format with:

```sh
git ls-files '*.h' '*.cpp' '*.c' | xargs clang-format-18 -i
```

[AGENTS.md](AGENTS.md) has the clang-tidy command, and the repository's
conventions for anyone, person or coding agent, making a change.

## Changes that move numbers

If your change can move a filter's estimates, covariance or innovation, say
so in the pull request, even when it is a fix. Say which filters and
configurations it affects and how you checked the size of the effect (a
test, a replay of recorded data, or one of the benchmarks under `examples/`).
The maintainers use this for the release notes' numerical-change entry.

## Interfaces

The interfaces covered by the [compatibility contract](RELEASING.md#the-compatibility-contract)
cannot change incompatibly outside a major release. If your change needs to,
open an issue first so that the maintainers can plan it; a compatible
alternative, such as a new overload or an opt-in parameter, is usually
possible. Code in `UKF::Detail` is not covered and can change freely.
