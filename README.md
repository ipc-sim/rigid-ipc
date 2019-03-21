# Fixing Collisions

A project for robustly resolving collisions with a guarantee of no interpenetrations and/or pass through.

* [Google Doc](https://docs.google.com/document/d/13MetSJoTTZ0ptT0SERbst1SgG-KbgK48hozhko6mJxc/edit?usp=sharing)
* [LaTeX Write-up](https://www.overleaf.com/6555952782nttqwfwgksjb)

## Compilation

To build the project, use the following commands from the root directory of the project.

```bash
mkdir build && cd build
cmake ..
make
```

### Dependencies

Fixing Collisions includes some option and non-optional dependencies.

#### NLopt

**ToDo:** add a build option to enable/disable NLopt
(**Build Option:** `-DENABLE_MOSEK=ON`)

NLopt is **automatically downloaded** through CMake.

NLopt is used for non-linear optimization of our constrained objective function. Currently, NLopt is required to make the project, but this may change in the near future.

#### Ipopt

**Build Option:** `-DENABLE_IPOPT=ON`

To install Ipopt via Homebrew, use the following commands.

```bash
brew tap udacity/CarND-MPC-Project https://github.com/udacity/CarND-MPC-Project
brew install ipopt --with-openblas
```

Ipopt is used for non-linear optimization of out constrained objective function. Ipopt is **not required** to make the project, but it is recommended in order to use the interior point optimization method.

#### OSQP

**Build Option:** `-DENABLE_OSQP=ON`

When OSQP is enabled for the first time, it is **downloaded automatically** through CMake.

OSQP is used for quadratic programming of the linearized constraints and the interior iterations of the nonlinear complementary problem. OSQP is not required to make the project, but it is recommended in order to use the linearized constraints.

#### MOSEK

**Build Option:** `-DENABLE_MOSEK=ON`

Currently, MOSEK must be installed manually (steps below).

MOSEK is used for quadratic programming of the linearized constraints and the interior iterations of the nonlinear complementary problem. MOSEK is **not required** to make the project, but it is recommended in order to use the linearized constraints.

*Currently, only MOSEK 7 works on macOS. MOSEK 8 gives the following error even after following all the installation instructions.*

```
dyld: Library not loaded: libmosek64.8.1.dylib
  Referenced from: <FIXING_COLLISIONS_DIR>/fixing-collisions/build/tests/unit_tests
  Reason: image not found
```

##### Installation Steps

1. Download MOSEK 7 from https://www.mosek.com/downloads/7.1.0.63/

2. Extract the files and follow MOSEK's installation directions found in `<MSKHOME>/mosek/quickstart.html` where `<MSKHOME>/` is the directory where MOSEK was extracted.

3. Get a personal academic license from https://www.mosek.com/products/academic-licenses/, and place it in the appropriate directory (`%USERPROFILE%\mosek\mosek.lic` on Windows and `$HOME/mosek/mosek.lic` on all UNIX like operating systems).

4. In order for CMake to find MOSEK you can either create a symbolic link:
```bash
ln -s <MSKHOME>/mosek /usr/local/mosek
```
or create a `MOSEK_DIR` environment variable:
```bash
export MOSEK_DIR="<MSKHOME>/mosek/7/tools/platform/osx64x86"
```

#### Catch2

Catch2 is **automatically downloaded** through CMake.

Catch2 is used for unit tests and is required to make the project.

#### libigl

Libigl must be **installed manually** (more information on downloading and compiling libigl can be found here https://libigl.github.io/).

Libigl is used for the viewer and for calling MOSEK's QP. Libigl is **required** to make the project.
