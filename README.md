# Rigid IPC

<img src="docs/imgs/teaser.png">

[![Build](https://github.com/ipc-sim/rigid-ipc/actions/workflows/continuous.yml/badge.svg)](https://github.com/ipc-sim/rigid-ipc/actions/workflows/continuous.yml)
[![License](https://img.shields.io/github/license/ipc-sim/rigid-ipc.svg?color=blue)](https://github.com/ipc-sim/rigid-ipc/blob/main/LICENSE)

<b>Robust, intersection-free, simulations of rigid bodies.</b>

This is the open-source reference implementation of the SIGGRAPH 2021 paper [Intersection-free Rigid Body Dynamics](https://ipc-sim.github.io/rigid-ipc/).

## Files

* `src/`: source code
* `cmake/` and `CMakeLists.txt`: CMake files
* `fixtures/`: input scripts to rerun all examples in our paper
* `meshes/`: input meshes used by the fixtures
* `tests/`: unit-tests
* `tools/`: Python and Bash scripts for generating and processing results
* `comparisons/`: files used in comparisons with other rigid body simulators
* `python/`: Python binding files
* `notebooks/`: Jupyter notebooks

## Build

To build the project, use the following commands from the root directory of the project:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```

### Dependencies

**All dependancies are downloaded through CMake** depending on the build options.
The following libraries are used in this project:

* [IPC Toolkit](https://github.com/ipc-sim/ipc-toolkit): common IPC functions
* [Eigen](https://eigen.tuxfamily.org/): linear algebra
* [libigl](https://github.com/libigl/libigl): basic geometry functions, predicates, and viewer
* [TBB](https://github.com/wjakob/tbb): parallelization
* [Tight Inclusion CCD](https://github.com/Continuous-Collision-Detection/Tight-Inclusion): correct (conservative) continuous collision detection between triangle meshes in 3D
* [spdlog](https://github.com/gabime/spdlog): logging information
* [filib](https://github.com/txstc55/filib): interval arithmetic
* [Niels Lohmann's JSON](https://github.com/nlohmann/json): parsing input JSON scenes
* [tinygltf](https://github.com/syoyo/tinygltf.git): exporting simulation animation to GLTF format
* [finite-diff](https://github.com/zfergus/finite-diff): finite difference comparisons
    * Only used by the unit tests and when `RIGID_IPC_WITH_DERIVATIVE_CHECK=ON`

#### Optional

* [Catch2](https://github.com/catchorg/Catch2.git): unit tests

## Scenes

We take as input a single JSON file that specifies the mesh and initial
conditions for each body. The `fixtures` directory contains example scenes.

## Python Bindings

We expose some functionality of Rigid IPC through Python. This is still in
development and lacks the ability to script many features available in the full
simulator.

To build the Python bindings use the `setup.py` script:
```sh
python setup.py install
```
