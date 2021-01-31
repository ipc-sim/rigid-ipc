# Rigid IPC

[![Build Status](https://travis-ci.com/geometryprocessing/fixing-collisions.svg?token=uiCkVsJXonpF1gk5xcRf&branch=master)](https://travis-ci.com/geometryprocessing/fixing-collisions)

**Robust, intersection free, simulations of rigid bodies in 2D and 3D.**

## Compilation

To build the project, use the following commands from the root directory of the project.

```bash
mkdir build
cd build
cmake ..
make -j4
```

### Dependencies

Most dependancies are downloaded through CMake depending on the build options.
The only exceptions to this are:

* [Boost](https://www.boost.org/): We currently use the interval arithmetic
library for interval root finding and the filesystem library for path
manipulation.

## Scenes

We take as input a single JSON file that specifies the mesh and initial
conditions for each body. The `fixtures` directory contains example scenes.
