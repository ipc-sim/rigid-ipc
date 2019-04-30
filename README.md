# Fixing Collisions

A project for robustly resolving collisions with a guarantee of no interpenetrations and/or pass through.

* [Google Doc](https://docs.google.com/document/d/13MetSJoTTZ0ptT0SERbst1SgG-KbgK48hozhko6mJxc/edit?usp=sharing)
* [LaTeX Write-up](https://www.overleaf.com/6555952782nttqwfwgksjb)

## TODO

* [x] Have one constraint per collision instead of per edge
* [ ] Line-search over each degree of freedom independently
    * [ ] Replace &gamma;, a scalar multiplier, with &Gamma;, a diagonal matrix multiplier.
* [ ] Update CCD to find all impacts where the time of impact is greater than or equal to t<sub>0</sub> not bounded by ≤ t<sub>1</sub>
    * The time of impact is recomputed every time we compute the constraints anyways
* [ ] Refactor OptimizationProblem
    * [ ] Compute **objective** value, gradient, and hessian in one function
    * [ ] Compute **constraint** value, gradient, and hessian in one function
* [ ] Add physics
    * [ ] Replace our objective with a more physically based one
    * [ ] Simple model (gravity, no friction, etc.)
    * [ ] Friction model
* [ ] Add static objects (displacement is fixed to zero)
    * [ ] Remove the variables from the optimization
    * OR
    * [ ] Set the displacements to zero every iteration (this might lead to intersections)
* [ ] Add rigid bodies by limiting the degrees of freedom to one rotation and position per body. Solve the optimization on these variables.
    * [ ] Implement translation from rigid motion to point-wise linear trajectories per time-step
* [ ] Add minimum separation distance
    * [ ] Trivially, change the lower bound of the barrier constraints

## Compilation

To build the project, use the following commands from the root directory of the project.

```bash
mkdir build && cd build
cmake ..
make
```

### Dependencies

Fixing Collisions includes some optional and non-optional dependencies.

#### NLopt

**ToDo:** add a build option to enable/disable NLopt
(**Build Option:** `-DENABLE_NLOPT=ON`)

NLopt is **automatically downloaded** through CMake.

NLopt is used for non-linear optimization of our constrained objective function. Currently, NLopt is required to make the project, but this may change in the near future.

#### Ipopt

**Build Option:** `-DENABLE_IPOPT=ON`

To install Ipopt via Homebrew, use the following commands.

```bash
brew tap udacity/CarND-MPC-Project https://github.com/udacity/CarND-MPC-Project
brew install ipopt --with-openblas
```

Ipopt is used for non-linear optimization of our constrained objective function. Ipopt is **not required** to make the project, but it is recommended in order to use the interior point optimization method.

#### OSQP

**Build Option:** `-DENABLE_OSQP=ON`

When OSQP is enabled for the first time, it is **downloaded automatically** through CMake.

OSQP is used for quadratic programming of the linearized constraints and the interior iterations of the nonlinear complementarity problem. OSQP is not required to make the project, but it is recommended in order to use the linearized constraints.

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

#### Catch2

**Build Option:** `-DBUILD_UNIT_TESTS=ON`

Catch2 is **automatically downloaded** through CMake.

Catch2 is used for unit tests. Catch2 is not required to make the project, but it is recommended in order to test the build.

#### libigl

Libigl must be **installed manually** (more information on downloading and compiling libigl can be found here https://libigl.github.io/).

Libigl is used for the viewer and for calling MOSEK's QP. Libigl is **required** to make the project.
