cmake_minimum_required(VERSION 3.1)
project(FixingCollisions)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_party/)

# libigl
option(LIBIGL_WITH_OPENGL               "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use IMGUI"          ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
set(SRCFILES
    src/collision_detection.cpp
    src/main.cpp
    src/state.cpp
    src/viewer.cpp
    src/viewer_menu.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
include_directories(${PROJECT_NAME}_bin include src)

add_library(json INTERFACE)
target_include_directories(json SYSTEM INTERFACE ${THIRD_PARTY_DIR}/json/include)
target_link_libraries(${PROJECT_NAME}_bin json)
