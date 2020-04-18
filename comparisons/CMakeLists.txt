cmake_minimum_required(VERSION 3.8)
################################################################################

option(FIXING_COLLISIONS_BUILD_BOX2D_COMPARISONS  "Build Box2D comparison"  ON)

if(FIXING_COLLISIONS_BUILD_BOX2D_COMPARISONS)
  add_subdirectory(Box2D)
endif()
