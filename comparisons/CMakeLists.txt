cmake_minimum_required(VERSION 3.8)
################################################################################

option(RIGID_IPC_WITH_BOX2D_COMPARISONS  "Build Box2D comparison"   ON)
option(RIGID_IPC_WITH_STIV_COMPARISONS   "Build Box2D comparison"  OFF)

if(RIGID_IPC_WITH_BOX2D_COMPARISONS)
  add_subdirectory(Box2D)
endif()

if(RIGID_IPC_WITH_STIV_COMPARISONS)
  add_subdirectory(STIV)
endif()
