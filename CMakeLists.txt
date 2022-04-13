# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(RIGID_IPC_TOPLEVEL_PROJECT OFF)
else()
    set(RIGID_IPC_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(RIGID_IPC_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Rigid IPC is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/RigidIPCOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/RigidIPCOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/RigidIPCOptions.cmake)
endif()

################################################################################

project(RigidIPC
        DESCRIPTION "Robust, intersection-free, simulations of rigid bodies."
        LANGUAGES CXX)

option(RIGID_IPC_WITH_UNIT_TESTS             "Build unit-tests"      ${RIGID_IPC_TOPLEVEL_PROJECT})
option(RIGID_IPC_WITH_OPENGL                 "Build GUI"             ${RIGID_IPC_TOPLEVEL_PROJECT})
option(RIGID_IPC_WITH_TOOLS                  "Build tools"           ${RIGID_IPC_TOPLEVEL_PROJECT})
option(RIGID_IPC_WITH_PROFILING              "Profile functions"                               OFF)
option(RIGID_IPC_WITH_COMPARISONS            "Build comparisons"                               OFF)
option(RIGID_IPC_WITH_SIMD                   "Enable SIMD"                                     OFF)
option(RIGID_IPC_WITH_PYTHON                 "Build Python bindings"                           OFF)
option(RIGID_IPC_WITH_DERIVATIVE_CHECK      "Check derivatives using finite differences"       OFF)

# Set default minimum C++ standard
if(RIGID_IPC_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/rigid_ipc/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# Color output
include(rigid_ipc_use_colors)

include(rigid_ipc_autogen)

mark_as_advanced(CLEAR CMAKE_BUILD_TYPE)

################################################################################
# External project settings
################################################################################

# libigl options
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_PREDICATES    "Use exact predicates"          ON)
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI ${RIGID_IPC_WITH_OPENGL} CACHE BOOL "Use OpenGL" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW       ${RIGID_IPC_WITH_OPENGL} CACHE BOOL "Use GLFW"   FORCE)
set(LIBIGL_WITH_OPENGL            ${RIGID_IPC_WITH_OPENGL} CACHE BOOL "Use IMGUI"  FORCE)
set(LIBIGL_WITH_PNG               ${RIGID_IPC_WITH_OPENGL} CACHE BOOL "Use PNG"    FORCE)

# PolySolve settings
option(POLYSOLVE_WITH_CHOLMOD          "Enable Cholmod library"            ON)
option(POLYSOLVE_WITH_UMFPACK          "Enable UmfPack library"            OFF)
option(POLYSOLVE_WITH_SUPERLU          "Enable SuperLU library"            OFF)
option(POLYSOLVE_WITH_MKL              "Enable MKL library"                OFF)
option(POLYSOLVE_WITH_PARDISO          "Enable Pardiso library"            OFF)
option(POLYSOLVE_WITH_HYPRE            "Enable hypre"                      OFF)
option(POLYSOLVE_WITH_AMGCL            "Use AMGCL"                         OFF)
option(POLYSOLVE_WITH_SPECTRA          "Enable computing spectrum"         OFF)

set(TIGHT_INCLUSION_WITH_NO_ZERO_TOI OFF CACHE BOOL "Enable refinement if CCD produces a zero ToI" FORCE)

################################################################################
# Rigid IPC Library
################################################################################

add_library(ipc_rigid
  src/autodiff/autodiff.cpp

  src/utils/sinc.cpp

  src/ccd/impact.cpp
  src/ccd/ccd.cpp
  src/ccd/linear/broad_phase.cpp
  src/ccd/piecewise_linear/time_of_impact.cpp
  src/interval/filib_rounding.cpp
  src/interval/interval_root_finder.cpp
  src/ccd/rigid/broad_phase.cpp
  src/ccd/rigid/rigid_body_hash_grid.cpp
  src/ccd/rigid/rigid_body_bvh.cpp
  src/ccd/rigid/time_of_impact.cpp
  src/ccd/rigid/rigid_trajectory_aabb.cpp
  src/ccd/redon/time_of_impact.cpp
  src/ccd/save_queries.cpp

  src/geometry/intersection.cpp

  src/io/serialize_json.cpp
  src/io/read_rb_scene.cpp
  src/io/read_obj.cpp
  src/io/write_obj.cpp
  src/io/write_gltf.cpp

  src/physics/mass.cpp
  src/utils/mesh_selector.cpp
  src/physics/rigid_body.cpp
  src/physics/rigid_body_assembler.cpp
  src/physics/rigid_body_problem.cpp

  src/barrier/barrier.cpp
  src/barrier/barrier_chorner.cpp

  src/opt/distance_barrier_constraint.cpp
  src/opt/collision_constraint.cpp
  src/opt/optimization_problem.cpp
  src/opt/optimization_results.cpp

  src/problems/barrier_problem.cpp

  src/solvers/newton_solver.cpp
  src/solvers/ipc_solver.cpp
  src/solvers/homotopy_solver.cpp
  src/solvers/solver_factory.cpp
  # src/solvers/line_search.cpp

  src/problems/rigid_body_collision_constraint.cpp
  src/problems/distance_barrier_rb_problem.cpp
  src/problems/split_distance_barrier_rb_problem.cpp
  src/problems/problem_factory.cpp

  src/time_stepper/sympletic_euler_time_stepper.cpp
  src/time_stepper/verlet_time_stepper.cpp
  src/time_stepper/exponential_euler_time_stepper.cpp
  src/time_stepper/dmv_time_stepper.cpp
  src/time_stepper/time_stepper_factory.cpp

  src/utils/tensor.cpp
  src/utils/eigen_ext.cpp
  src/utils/regular_2d_grid.cpp
  src/utils/get_rss.cpp

  src/SimState.cpp
  src/logger.cpp
  src/profiler.cpp
)
target_include_directories(ipc_rigid PUBLIC src)
add_library(ipc::rigid ALIAS ipc_rigid)

# Autogenerated files ----------------------------------------------------------
# project_autogen(ipc_rigid time_of_impact_coeff.py time_of_impact_coeff.tpp auto_time_of_impact_coeff .ipp)

################################################################################
# Optional Definitions
################################################################################

# Directory paths
set(RIGID_IPC_MESH_DIR ${CMAKE_CURRENT_SOURCE_DIR}/meshes)
target_compile_definitions(ipc_rigid PUBLIC RIGID_IPC_MESHES_DIR=\"${RIGID_IPC_MESH_DIR}\")

if(RIGID_IPC_WITH_DERIVATIVE_CHECK)
  target_compile_definitions(ipc_rigid PUBLIC RIGID_IPC_WITH_DERIVATIVE_CHECK)
endif()

if(RIGID_IPC_WITH_PROFILING)
  message(STATUS "Profiling Enabled")
  target_compile_definitions(ipc_rigid PUBLIC RIGID_IPC_PROFILE_FUNCTIONS)
endif()

# For MSVC, do not use the min and max macros.
target_compile_definitions(ipc_rigid PUBLIC NOMINMAX)

################################################################################
# Required Libraries
################################################################################

# Extra warnings
include(rigid_ipc_warnings)
target_link_libraries(ipc_rigid PRIVATE ipc::rigid::warnings)

# libigl
include(eigen)
include(libigl)
target_link_libraries(ipc_rigid PUBLIC
  Eigen3::Eigen
  igl::core
  igl::predicates
)

# JSON Parser
include(json)
target_link_libraries(ipc_rigid PUBLIC nlohmann::json)

# Logger
include(spdlog)
target_link_libraries(ipc_rigid PUBLIC spdlog::spdlog)

# Finite Differences
include(finite_diff)
target_link_libraries(ipc_rigid PUBLIC finitediff::finitediff)

# Boost
include(boost)
target_link_libraries(ipc_rigid PUBLIC Boost::boost)

# TBB
include(tbb)
target_link_libraries(ipc_rigid PUBLIC TBB::tbb)

# MSCCD
include(tight_inclusion)
target_link_libraries(ipc_rigid PUBLIC tight_inclusion::tight_inclusion)

# PolySolve for wrapping linear solvers
include(polysolve)
target_link_libraries(ipc_rigid PUBLIC PolyFEM::polysolve)

# IPC Toolkit
include(ipc_toolkit)
target_link_libraries(ipc_rigid PUBLIC ipc::toolkit)

# filib
include(filib)
target_link_libraries(ipc_rigid PUBLIC filib::filib)

# SimpleBVH
include(simple_bvh)
target_link_libraries(ipc_rigid PUBLIC simple_bvh::simple_bvh)

# tinygltf
include(tinygltf)
target_link_libraries(ipc_rigid PUBLIC tinygltf::tinygltf)

# Filesystem
include(filesystem)
target_link_libraries(ipc_rigid PUBLIC ghc::filesystem)

################################################################################
# Compiler options
################################################################################

## SIMD support
if(RIGID_IPC_WITH_SIMD)
  # Figure out SSE level support
  message(STATUS "Seaching for SSE...")
  find_package(SSE)
  # Apple broke AVX instructions with clang:
  # https://gitlab.com/libeigen/eigen/-/commit/71aa53dd6dfdc497324d9e87f59c4ba820191856
  if(NOT APPLE)
    # Figure out AVX level support
    message(STATUS "Searching for AVX...")
    find_package(AVX)
    # Figure out FMA level support
    message(STATUS "Searching for FMA...")
    find_package(FMA)
  else()
    message(WARNING "Disabling AVX support: clang compiler shipped with XCode 11.[012] generates broken assembly with -macosx-version-min=10.15 and AVX enabled.")
  endif()
  # Add SSE, AVX, and FMA flags to compiler flags
  string(REPLACE " " ";" SIMD_FLAGS "${SSE_FLAGS} ${AVX_FLAGS} ${FMA_FLAGS}")
  target_compile_options(ipc_rigid PUBLIC ${SIMD_FLAGS})
endif()

# Use C++17
target_compile_features(ipc_rigid PUBLIC cxx_std_17)

################################################################################
# Simulations Executable
################################################################################

add_executable(rigid_ipc_sim src/main.cpp)
if(RIGID_IPC_WITH_OPENGL)
    target_sources(rigid_ipc_sim PUBLIC
      src/viewer/imgui_ext.cpp
      src/viewer/igl_viewer_ext.cpp
      src/viewer/UISimState.cpp
      src/viewer/UIMenu.cpp
    )
endif()

target_link_libraries(rigid_ipc_sim PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(rigid_ipc_sim PUBLIC CLI11::CLI11)

if(RIGID_IPC_WITH_OPENGL)
    target_compile_definitions(rigid_ipc_sim PUBLIC RIGID_IPC_WITH_OPENGL)

    include(libigl)
    target_link_libraries(rigid_ipc_sim PUBLIC
        igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::png)

    # Charlie Tangora's gif-h library
    include(gif_h)
    target_link_libraries(rigid_ipc_sim PRIVATE gif_h::gif_h)
endif()

################################################################################
# Tools
################################################################################
if(RIGID_IPC_WITH_TOOLS)
  add_subdirectory(src/tools)
endif()

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(RIGID_IPC_WITH_UNIT_TESTS)
    include(CTest)
    enable_testing()

    # Include Catch2 and provide function `catch_discover_tests` to register tests.
    include(catch2)
    FetchContent_GetProperties(catch2)
    include("${catch2_SOURCE_DIR}/contrib/Catch.cmake")

    add_subdirectory(tests)
endif()

################################################################################
# Comparisons to other methods
################################################################################
if(RIGID_IPC_WITH_COMPARISONS)
  add_subdirectory(comparisons)
endif()

################################################################################

if(SKBUILD)
  # Scikit-Build does not add your site-packages to the search path
  # automatically, so we need to add it _or_ the pybind11 specific directory
  # here.
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "import pybind11; print(pybind11.get_cmake_dir())"
    OUTPUT_VARIABLE _tmp_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
  list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
endif()

if (UNIX AND NOT APPLE)
    add_compile_options(-std=c++1y)
    add_compile_options(-pthread)
endif()

include(pybind11)

# find_package(pybind11 CONFIG REQUIRED)
# find_package(PhoXi REQUIRED CONFIG PATHS "$ENV{PHOXI_CONTROL_PATH}")

pybind11_add_module(_${PROJECT_NAME} MODULE src/bindings.cpp)
target_link_libraries(_${PROJECT_NAME} PUBLIC ipc::rigid)

# target_include_directories(_${PROJECT_NAME} PUBLIC ${PHOXI_INCLUDE_DIRS})
# target_link_libraries(_${PROJECT_NAME}
#   PRIVATE
#   ${PHOXI_LIBRARY}
#   $<$<PLATFORM_ID:Linux>:rt>
# )
target_compile_definitions(_${PROJECT_NAME} PRIVATE VERSION_INFO=${PROJECT_VERSION})

set(lib_path "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${lib_path}" is_system)
if ("${is_system}" STREQUAL "-1")
    set_target_properties(_${PROJECT_NAME} PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_RPATH "${lib_path}")
endif()

install(TARGETS _${PROJECT_NAME}
    DESTINATION .
)