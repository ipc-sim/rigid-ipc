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
        GIT_TAG        v2.5.0
    )
endfunction()

## nlopt
function(download_nlopt)
    custom_download_project(nlopt
        GIT_REPOSITORY https://github.com/stevengj/nlopt.git
        GIT_TAG        463abb855d166503a4f2890379647a54a26ca678
    )
endfunction()

## OSQP
function(download_osqp)
    custom_download_project(osqp
        GIT_REPOSITORY https://github.com/oxfordcontrol/osqp.git
        GIT_TAG        8949e678122d6949139644d2a95985765527535f
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

function(download_libigl)
    custom_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        aea868bd1fc64f71afecd2c51e51507a99d8e3e5
    )
endfunction()

function(download_json)
    custom_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json
        GIT_TAG v3.7.0
    )
endfunction()
