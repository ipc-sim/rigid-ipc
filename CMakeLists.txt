cmake_minimum_required(VERSION 3.0)
project(FixingCollisions)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# include
include_directories(include src)

# target
# TODO: Add a real executable target here (e.g. the GUI)
add_executable(fixing-collisions src/collision_detection.cpp)

# external libraries
# add_subdirectory(libs/...)
target_link_libraries(fixing-collisions Eigen3::Eigen)
