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
        GIT_TAG        v2.11.3
    )
endfunction()

function(fixing_collisions_download_fmt)
     fixing_collisions_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        6.2.0
    )
endfunction()

function(fixing_collisions_download_spdlog)
     fixing_collisions_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.5.0
    )
endfunction()

function(fixing_collisions_download_libigl)
     fixing_collisions_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        v2.2.0
    )
endfunction()

function(fixing_collisions_download_json)
     fixing_collisions_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.7.3
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
        GIT_TAG        dac4e629defb0edef050e56050eefb00574ecbaa
    )
endfunction()

function(fixing_collisions_download_tbb)
   fixing_collisions_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        20357d83871e4cb93b2c724fe0c337cd999fd14f
  )
endfunction()

# Wrapper for various CCD codes
# function(fixing_collisions_download_ccd_wrapper)
#   fixing_collisions_download_project(ccd-wrapper
#     GIT_REPOSITORY https://github.com/zfergus/ccd-wrapper.git
#     GIT_TAG        3e1eada7ae6e04fdbe4572417fc05c0b4f698cd6
#   )
# endfunction()

# Etienne Vouga's CTCD Library
function(fixing_collisions_download_evctcd)
  fixing_collisions_download_project(EVCTCD
    GIT_REPOSITORY https://github.com/evouga/collisiondetection.git
    GIT_TAG        e5fe5c9767207df5047e375fb20180a665ae186f
  )
endfunction()

# Rational CCD (rational version of Brochu et al. [2012])
function(fixing_collisions_download_rational_ccd)
  fixing_collisions_download_project(rational_ccd
    GIT_REPOSITORY https://github.com/teseoch/Exact-CDD.git
    GIT_TAG        63fdb92b917c8d7349c9403656f1497df6fa6ce5
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

# Polysolve for Linear Solver Wrapper
function(fixing_collisions_download_polysolve)
  fixing_collisions_download_project(polysolve
    GIT_REPOSITORY https://github.com/polyfem/polysolve.git
    GIT_TAG        23ccf9487772ab1b9554a78d8808179d1b3286a0
  )
endfunction()

# IPC Toolkit
function(fixing_collisions_download_ipc_toolkit)
  fixing_collisions_download_project(ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG        521366a6c6c756e4f92c55d6514654376d9b494e
  )
endfunction()
