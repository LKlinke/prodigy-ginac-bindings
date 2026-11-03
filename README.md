# prodigy-ppl-ginac: A library for manipulating discrete distributions through generating functions

GiNaC-based Python bindings used by [Prodigy](https://github.com/LKlinke/Prodigy) for its fast `ginac` engine.

The main features are
- python bindings
- a purely string-based interface, no need to mess around with data structures for algebraic expression and the like
- support for parameters (can also be used as meta symbols in second-order PGFs)

## Installation

```
pip install prodigy-ppl-ginac
```

Prebuilt wheels contain GiNaC, CLN and GMP, so they need neither a compiler nor a system GiNaC.
They exist for CPython 3.11 to 3.14 on
- Linux (manylinux_2_28) x86_64 and aarch64,
- macOS 11 or later, x86_64 and arm64.

Windows is not supported.

## Usage

The distribution is called `prodigy-ppl-ginac`, the module `prodigy_ginac`:

```
import prodigy_ginac
dist = prodigy_ginac.geometric("x", "1/2")
e = dist.E("x")  # e becomes 1, the expected value of the 1/2-geometric distribution
```

Before it was published on PyPI, the module was installed from Git as `pygin` (version 1.0); the package `pygin` on PyPI
is an unrelated project.

## Building from source

On other platforms pip builds the sdist, which needs a C++17 compiler, CMake and a system GiNaC with pkg-config,
e.g. `brew install ginac pkgconf` on macOS or `apt install libginac-dev pkg-config` on Debian and Ubuntu.
In a clone of this repository

```
pip install .
```

builds and installs the module.

`scripts/build-deps.sh PREFIX` builds the pinned versions of GMP, CLN and GiNaC as static libraries into `PREFIX`,
as for the wheels; with `CMAKE_PREFIX_PATH=PREFIX` the build uses them instead of the system GiNaC.

### Wheels and releases

[cibuildwheel](https://cibuildwheel.pypa.io) builds the wheels (configuration in `pyproject.toml`, workflow
`.github/workflows/wheels.yml`). Pushing a tag `v<version>` (e.g. `v0.0.1`, or `v0.1.0rc1` for a pre-release) publishes the wheels and the sdist with
this version on PyPI (`.github/workflows/release.yml`). The version is in no file; without a tag, builds get a
development version such as `0.1.devN`.

### C++ tests and demo

The C++ tests need the googletest submodule:

```
git submodule update --init
cmake -B build && cmake --build build
./build/test/prodigy_ginac_test
```

The CMake options `PRODIGY_GINAC_BUILD_TESTS` and `PRODIGY_GINAC_BUILD_DEMO` (both on by default, off for the wheels)
switch the tests and the demo executable.

## License

The source code in this repository is licensed under the Apache License 2.0 (`LICENSE`).

The wheels on PyPI also contain statically linked copies of
[GiNaC](https://www.ginac.de) and [CLN](https://www.ginac.de/CLN/) (GPL-2.0-or-later) and
[GMP](https://gmplib.org) (LGPL-3.0-or-later or GPL-2.0-or-later). A wheel as a whole is therefore distributed under
the terms of the GNU General Public License, version 3 or later, the GPL version compatible with the Apache License 2.0.
The license texts of these libraries are in the wheel's `.dist-info/licenses/third_party` directory; their source
code is available from the links above, in the versions pinned in `scripts/build-deps.sh`.
