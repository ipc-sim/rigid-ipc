include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FIXING_COLLISIONS_EXTRA_OPTIONS "")
endif()

function(fixing_collisions_download_project name)
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

function(fixing_collisions_download_catch2)
     fixing_collisions_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.5.0
    )
endfunction()

function(fixing_collisions_download_fmt)
     fixing_collisions_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        5.3.0
    )
endfunction()

function(fixing_collisions_download_spdlog)
     fixing_collisions_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.3.1
    )
endfunction()

function(fixing_collisions_download_libigl)
     fixing_collisions_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        aea868bd1fc64f71afecd2c51e51507a99d8e3e5
    )
endfunction()

function(fixing_collisions_download_json)
     fixing_collisions_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json
        GIT_TAG v3.7.0
    )
endfunction()

function(fixing_collisions_download_cli11)
     fixing_collisions_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        8ecce8fd2c49f64c80e5757cb12d2fd1fa62f242
    )
endfunction()

function(fixing_collisions_download_finite_diff)
     fixing_collisions_download_project(finite-diff
        GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
        GIT_TAG        01e16fb39ecf6fe6367b932448028817becfb633
    )
endfunction()

function(fixing_collisions_download_tbb)
   fixing_collisions_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        20357d83871e4cb93b2c724fe0c337cd999fd14f
  )
endfunction()

# Wrapper for various CCD codes
function(fixing_collisions_download_ccd_wrapper)
  fixing_collisions_download_project(ccd-wrapper
    GIT_REPOSITORY https://github.com/zfergus/ccd-wrapper.git
    GIT_TAG        4e86dea1a37b5899c01a9bd2cbf2744769885097
  )
endfunction()

# xml for MJCF format
function(fixing_collisions_download_tinyxml2)
   fixing_collisions_download_project(tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG        bfbcc0333d1a24ec8d9e10d14116d00dbdedf043
  )
endfunction()

# GIF writer
function(fixing_collisions_download_gif_h)
   fixing_collisions_download_project(gif-h
    GIT_REPOSITORY https://github.com/charlietangora/gif-h.git
    GIT_TAG        606168840bff429a894dd470e36fb9956cbb07c0
  )
endfunction()
