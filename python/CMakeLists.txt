cmake_minimum_required(VERSION 3.8)

# pybind11
include(pybind11)

# Create Python module for Rigid IPC
pybind11_add_module(rigidipc src/bindings.cpp)
target_link_libraries(rigidipc PRIVATE ipc_rigid)
# set_target_properties(rigidipc PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
