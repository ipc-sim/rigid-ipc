# Fixing Collisions

[![Build Status](https://travis-ci.com/geometryprocessing/fixing-collisions.svg?token=uiCkVsJXonpF1gk5xcRf&branch=master)](https://travis-ci.com/geometryprocessing/fixing-collisions)

**Robust, intersection free, simulations of rigid bodies in 2D and 3D.**

* [LaTeX Write-up](https://www.overleaf.com/6555952782nttqwfwgksjb)
<!-- * [Google Doc](https://docs.google.com/document/d/13MetSJoTTZ0ptT0SERbst1SgG-KbgK48hozhko6mJxc/edit?usp=sharing) -->

## Compilation

To build the project, use the following commands from the root directory of the project.

```bash
mkdir build
cd build
cmake ..
make
```

### Dependencies

Most dependancies are downloaded through CMake depending on the build options.
The only exceptions to this are:

* [Boost](https://www.boost.org/): We currently use the interval arithmetic
library for interval root finding and the filesystem library for path
manipulation.
* Python: We use SymPy to generate some code automatically. See
`requirements.txt` for a list of required packages.

<!--
#### MOSEK

**Build Option:** `-DENABLE_MOSEK=ON`

Currently, MOSEK must be installed manually (steps below).

MOSEK is used for quadratic programming of the linearized constraints and the interior iterations of the nonlinear complementarity problem. MOSEK is **not required** to make the project, but it is recommended in order to use the linearized constraints.

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
-->

## Scenes

We take as input a single JSON file that specifies the mesh and initial
conditions for each body. The `fixtures` directory contains example scenes.

## Algorithm

```
x_0^{t+1} = x^{t} + h * v^{t} + h^2 * g
v_0^{t+1} = (x_0^{t+1} - x^{t}) / h
if (collision_between(x^{t}, x_0^{t+1})) {
    original_collisions = CCD(x^{t}, x_0^{t+1})
    x^{t+1} = barrier_newton_solve(x^{t}, x_0^{t+1})
    if (using sequential impulses for restitution) {
         // Apply the restitution model in the above write-up
        v^{t+1} = solve_velocities(
            x^{t}, x_0^{t+1},
            v^{t}, v_0^{t+1},
            original_collisions,
            coefficient_of_restitution)
    } else {
        v^{t+1} = x^{t+1} - x^{t} / h
    }
} else {
    x^{t+1} = x_0^{t+1}
    v^{t+1} = v_0^{t+1}
}
```
