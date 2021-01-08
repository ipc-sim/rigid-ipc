# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(FixingCollisionsDownloadExternal)

################################################################################
# Required libraries
################################################################################

# MOSEK library
if(FIXING_COLLISIONS_ENABLE_MOSEK)
    if(NOT TARGET mosek)
        # fixing_collisions_download_mosek()
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
endif()

if(NOT TARGET fmt::fmt)
    fixing_collisions_download_fmt()
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
    fixing_collisions_download_spdlog()
    add_library(spdlog INTERFACE)
    add_library(spdlog::spdlog ALIAS spdlog)
    target_include_directories(spdlog SYSTEM INTERFACE ${FIXING_COLLISIONS_EXTERNAL}/spdlog/include)
    target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
    target_link_libraries(spdlog INTERFACE fmt::fmt)
endif()

# libigl
if(NOT TARGET igl::core)
    fixing_collisions_download_libigl()
    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${FIXING_COLLISIONS_EXTERNAL}/libigl/cmake")
    include(libigl)
endif()

# json
if(NOT TARGET nlohmann_json::nlohmann_json)
    fixing_collisions_download_json()
    option(JSON_BuildTests "" OFF)
    option(JSON_MultipleHeaders "" ON)
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/json json)
endif()

if(NOT TARGET CLI11::CLI11)
    fixing_collisions_download_cli11()
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/cli11)
endif()

# finite-diff
if(NOT TARGET FiniteDiff::FiniteDiff)
  fixing_collisions_download_finite_diff()
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/finite-diff EXCLUDE_FROM_ALL)
  add_library(FiniteDiff::FiniteDiff ALIAS FiniteDiff)
endif()

# TBB
if(NOT TARGET TBB::tbb)
  fixing_collisions_download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  add_library(TBB::tbb ALIAS tbb_static)
endif()

if(NOT TARGET CCDWrapper)
  fixing_collisions_download_ccd_wrapper()
  set(CCD_WRAPPER_WITH_UNIT_TESTS OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_BENCHMARK OFF CACHE BOOL " " FORCE)
  # methods
  # This is used for linearized CCD
  set(CCD_WRAPPER_WITH_FPRF ON CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_MSRF OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_RP OFF CACHE BOOL " " FORCE)
  # This is used for exact intersection check
  set(CCD_WRAPPER_WITH_RRP ON CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_FIXEDRP OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_BSC OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_TIGHT_CCD OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_INTERVAL OFF CACHE BOOL " " FORCE)
  set(CCD_WRAPPER_WITH_TIGHT_INCLUSION ON CACHE BOOL " " FORCE)
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/ccd-wrapper)
endif()

if(NOT TARGET tinyxml2::tinyxml2)
  fixing_collisions_download_tinyxml2()
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/tinyxml2)
  add_library(tinyxml2::tinyxml2 ALIAS tinyxml2)
endif()

if(NOT TARGET gif_h::gif_h)
  fixing_collisions_download_gif_h()
  add_library(gif_h INTERFACE)
  target_include_directories(gif_h INTERFACE ${FIXING_COLLISIONS_EXTERNAL}/gif-h)
  add_library(gif_h::gif_h ALIAS gif_h)
endif()

# Polysolve for Linear Solver Wrapper
if(NOT TARGET polysolve)
  fixing_collisions_download_polysolve()
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/polysolve)
endif()

# IPC Toolkit
if(NOT TARGET IPCToolkit)
  fixing_collisions_download_ipc_toolkit()
  set(IPC_TOOLKIT_BUILD_UNIT_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/ipc-toolkit)
endif()

# filib
if(NOT TARGET filib)
  fixing_collisions_download_filib()
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/filib)
endif()

# SimpleBVH
if(NOT TARGET SimpleBVH::BVH_lib)
  fixing_collisions_download_simple_bvh()
  add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/SimpleBVH)
  add_library(SimpleBVH::BVH_lib ALIAS BVH_lib)
endif()
