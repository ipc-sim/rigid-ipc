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

## nlopt
function(download_nlopt)
    custom_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        master
    )
endfunction()

## OSQP
function(download_osqp)
    custom_download_project(osqp
        GIT_REPOSITORY https://github.com/oxfordcontrol/osqp.git
        GIT_TAG        master
    )
endfunction()

function(download_fmt)
    custom_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        5.3.0
    )
endfunction()

function(download_spdlog)
    custom_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.3.1
    )
endfunction()

## libigl
## GIT_TAG dev
## GIT_HASH aea868bd1fc64f71afecd2c51e51507a99d8e3e5
