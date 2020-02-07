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

# NLopt library
if(FIXING_COLLISIONS_ENABLE_NLOPT)
    if(NOT TARGET nlopt)
        fixing_collisions_download_nlopt()
        find_package(NLopt QUIET)
        if(NLOPT_FOUND)
            message(STATUS "Including NLOPT")
            add_library(nlopt_nlopt INTERFACE)
            target_link_libraries(nlopt_nlopt INTERFACE ${NLOPT_LIBRARIES})
            add_library(nlopt::nlopt ALIAS nlopt_nlopt)
        else()
            message(WARNING "NLopt not found!")
            add_library(nlopt::nlopt INTERFACE IMPORTED)
        endif()
    endif()
endif()

# OSQP library
if(FIXING_COLLISIONS_ENABLE_OSQP)
    if(NOT TARGET osqp::osqp)
        fixing_collisions_download_osqp()
        # Make sure the right types are used
        set(DFLOAT OFF CACHE BOOL "Use float numbers instead of doubles"   FORCE)
        set(DLONG  OFF CACHE BOOL "Use long integers (64bit) for indexing" FORCE)
        add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/osqp EXCLUDE_FROM_ALL)
        if(UNIX AND NOT APPLE)
            set_target_properties(osqpstatic PROPERTIES INTERFACE_LINK_LIBRARIES ${CMAKE_DL_LIBS})
        endif()
        add_library(osqp::osqp ALIAS osqpstatic)
    endif()
endif()

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
