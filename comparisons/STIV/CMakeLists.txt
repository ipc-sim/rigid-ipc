cmake_minimum_required(VERSION 3.8)

################################################################################

option(STIV_WITH_MOSEK "Enable Mosek Library" OFF)
mark_as_advanced(STIV_WITH_MOSEK)

################################################################################

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

################################################################################

add_library(STIV

  src/opt/volume_constraint.cpp

  src/solvers/fischer_newton.cpp
  src/solvers/lcp_solver.cpp
  src/solvers/ncp_solver.cpp
  src/solvers/line_search.cpp

  src/problems/volume_rb_problem.cpp
)

target_include_directories(STIV PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# project_autogen(STIV collision_volume.py collision_volume.tpp auto_collision_volume .ipp)

# MOSEK QP solver
if(STIV_WITH_MOSEK)
  if(NOT TARGET mosek)
    find_package(MOSEK QUIET)
    if(MOSEK_FOUND)
      message(STATUS "Including MOSEK")
      # Make sure libigl uses mosek
      set(LIBIGL_WITH_MOSEK ON CACHE BOOL "Use MOSEK" FORCE)
      # Create a library for mosek
      add_library(mosek_mosek INTERFACE)
      target_link_libraries(mosek_mosek INTERFACE ${MOSEK_LIBRARIES})
      target_include_directories(mosek_mosek SYSTEM INTERFACE ${MOSEK_INCLUDE_DIRS})
      target_compile_definitions(mosek_mosek INTERFACE -DHAS_MOSEK)
      add_library(mosek::mosek ALIAS mosek_mosek)
    else()
      message(WARNING "MOSEK not found!")
      add_library(mosek::mosek INTERFACE IMPORTED)
    endif()
  endif()

  target_compile_definitions(${PROJECT_NAME} PUBLIC -DBUILD_WITH_MOSEK)
  target_link_libraries(${PROJECT_NAME} PUBLIC mosek::mosek)
endif()
