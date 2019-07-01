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
if(NOT TARGET nlopt)
    download_nlopt()
    find_package(NLopt QUIET)
    if(NLOPT_FOUND)
        message(STATUS "Including NLOPT")
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
    if(NOT TARGET osqp::osqp)
        download_osqp()
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
if(ENABLE_MOSEK)
    if(NOT TARGET mosek)
        # download_mosek()
        find_package(mosek QUIET)
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
    download_fmt()
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
    download_spdlog()
    add_library(spdlog INTERFACE)
    add_library(spdlog::spdlog ALIAS spdlog)
    target_include_directories(spdlog SYSTEM INTERFACE ${FIXING_COLLISIONS_EXTERNAL}/spdlog/include)
    target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
    target_link_libraries(spdlog INTERFACE fmt::fmt)
endif()

# libigl
if(NOT TARGET igl::core)
    download_libigl()

    # Import libigl targets
    list(APPEND CMAKE_MODULE_PATH "${FIXING_COLLISIONS_EXTERNAL}/libigl/cmake")
    include(libigl)
endif()

# json
if(NOT TARGET nlohmann_json::nlohmann_json)
    download_json()
    option(JSON_BuildTests "" OFF)
    option(JSON_MultipleHeaders "" ON)
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/json json)
endif()

# LBFGS++
# if(NOT TARGET lbfgspp::lbfgspp)
#     download_lbfgspp()
#     add_library(lbfgspp INTERFACE)
#     target_include_directories(lbfgspp SYSTEM INTERFACE ${FIXING_COLLISIONS_EXTERNAL}/lbfgspp/include)
#     add_library(lbfgspp::lbfgspp ALIAS lbfgspp)
# endif()
