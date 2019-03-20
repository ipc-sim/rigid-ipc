# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.

### Configuration
set(FIXING_COLLISIONS_ROOT     "${CMAKE_CURRENT_LIST_DIR}/..")
set(FIXING_COLLISIONS_EXTERNAL "${FIXING_COLLISIONS_ROOT}/3rd_party")

# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(FixingCollisionsDownloadExternal)

################################################################################
# Required libraries
################################################################################

# NLopt library
if(NOT TARGET nlopt)
    download_nlopt()
    find_package(NLopt QUIET)
    if(NLOPT_FOUND)
        message("Including NLOPT")
        add_library(nlopt_nlopt INTERFACE)
        target_link_libraries(nlopt_nlopt INTERFACE ${NLOPT_LIBRARIES})
        target_compile_definitions(nlopt_nlopt INTERFACE -DHAS_NLOPT)
        add_library(nlopt::nlopt ALIAS nlopt_nlopt)
    else()
        message(WARNING "NLopt not found!")
        add_library(nlopt::nlopt INTERFACE IMPORTED)
    endif()
endif()

# OSQP library
if(ENABLE_OSQP)
    if(NOT TARGET osqp)
        download_osqp()
        find_package(osqp QUIET)
        if(OSQP_FOUND)
            message("Including OSQP")
            # Make sure the right types are used
            set(DFLOAT OFF CACHE BOOL "Use float numbers instead of doubles"   FORCE)
            set(DLONG  OFF CACHE BOOL "Use long integers (64bit) for indexing" FORCE)
            # Create a library for OSQP
            add_library(osqp_osqp INTERFACE)
            target_link_libraries(osqp_osqp INTERFACE ${OSQP_LIBRARIES})
            target_compile_definitions(osqp_osqp INTERFACE -DHAS_OSQP)
            add_library(osqp::osqp ALIAS osqp_osqp)
        else()
            message(WARNING "OSQP not found!")
            add_library(osqp::osqp INTERFACE IMPORTED)
        endif()
    endif()
endif()

# Mosek library
if(ENABLE_MOSEK)
    if(NOT TARGET mosek)
        # download_mosek()
        find_package(mosek QUIET)
        if(MOSEK_FOUND)
            message("Including Mosek")
            # Make sure libigl uses mosek
            set(LIBIGL_WITH_MOSEK ON CACHE BOOL "Use MOSEK" FORCE)
            # Create a library for mosek
            add_library(mosek_mosek INTERFACE)
            target_link_libraries(mosek_mosek INTERFACE ${MOSEK_LIBRARIES})
            target_include_directories(mosek_mosek INTERFACE ${MOSEK_INCLUDE_DIRS})
            target_compile_definitions(mosek_mosek INTERFACE -DHAS_MOSEK)
            add_library(mosek::mosek ALIAS mosek_mosek)
        else()
            message(WARNING "Mosek not found!")
            add_library(mosek::mosek INTERFACE IMPORTED)
        endif()
    endif()
endif()
