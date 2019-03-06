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
        add_library(nlopt_nlopt INTERFACE)
        target_link_libraries(nlopt_nlopt INTERFACE ${NLOPT_LIBRARIES})
        target_compile_definitions(nlopt_nlopt INTERFACE -DHAS_NLOPT)
        add_library(nlopt::nlopt ALIAS nlopt_nlopt)
    else()
        message(STATUS "NLopt not found!")
        add_library(nlopt::nlopt INTERFACE IMPORTED)
    endif()
endif()
