# Rigid IPC

<!-- [![Build status](https://github.com/ipc-sim/IPC/workflows/Build/badge.svg?event=push)](https://github.com/ipc-sim/IPC/actions?query=workflow%3ABuild+branch%3Amaster+event%3Apush) -->
<!-- [![License](https://img.shields.io/github/license/ipc-sim/IPC.svg?color=blue)](https://github.com/ipc-sim/IPC/blob/master/LICENSE) -->

<img src="docs/imgs/teaser.png">

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

Most dependancies are downloaded through CMake depending on the build options.
The only exceptions to this are:

* [Boost](https://www.boost.org/): We currently use the interval arithmetic
library for interval root finding

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
