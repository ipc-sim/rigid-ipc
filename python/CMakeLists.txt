cmake_minimum_required(VERSION 3.8)

###############################################################################
# Download pybind11
function(rigid_ipc_download_pybind11)
  rigid_ipc_download_project(pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.6.2
  )
endfunction()

if(NOT TARGET pybind11::module)
  rigid_ipc_download_pybind11()
  add_subdirectory(
    ${RIGID_IPC_EXTERNAL}/pybind11
    ${CMAKE_CURRENT_BINARY_DIR}/pybind11)
endif()

###############################################################################
# Create Python module for Rigid IPC

pybind11_add_module(rigidipc src/bindings.cpp)
target_link_libraries(rigidipc PRIVATE ${PROJECT_NAME})
# set_target_properties(rigidipc PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
