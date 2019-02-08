include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "")
endif()

function(custom_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FIXING_COLLISIONS_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FIXING_COLLISIONS_EXTERNAL}/.cache/${name}
        QUIET
        ${FIXING_COLLISIONS_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################


## Catch2
function(download_catch2)
    custom_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.4.2
    )
endfunction()
