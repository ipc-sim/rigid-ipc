include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(RIGID_IPC_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(RIGID_IPC_EXTRA_OPTIONS "")
endif()

function(rigid_ipc_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${RIGID_IPC_EXTERNAL}/${name}
        DOWNLOAD_DIR ${RIGID_IPC_EXTERNAL}/.cache/${name}
        QUIET
        ${RIGID_IPC_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

function(rigid_ipc_download_fmt)
     rigid_ipc_download_project(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        6.2.0
    )
endfunction()

function(rigid_ipc_download_spdlog)
     rigid_ipc_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.5.0
    )
endfunction()

function(rigid_ipc_download_libigl)
     rigid_ipc_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        efee81b7dbc81ec87adaca1197b47f4faab961d3
    )
endfunction()

function(rigid_ipc_download_json)
     rigid_ipc_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.7.3
    )
endfunction()

function(rigid_ipc_download_cli11)
     rigid_ipc_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        8ecce8fd2c49f64c80e5757cb12d2fd1fa62f242
    )
endfunction()

function(rigid_ipc_download_finite_diff)
     rigid_ipc_download_project(finite-diff
        GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
        GIT_TAG        dac4e629defb0edef050e56050eefb00574ecbaa
    )
endfunction()

function(rigid_ipc_download_tbb)
   rigid_ipc_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        141b0e310e1fb552bdca887542c9c1a8544d6503
  )
endfunction()

# Tight Inclusion for MSCCD
function(rigid_ipc_download_tight_inclusion)
   rigid_ipc_download_project(Tight-Inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG        8730ca944b71b8bed6f1aecb0fca5bd152f95a58
  )
endfunction()

# xml for MJCF format
function(rigid_ipc_download_tinyxml2)
   rigid_ipc_download_project(tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG        bfbcc0333d1a24ec8d9e10d14116d00dbdedf043
  )
endfunction()

# GIF writer
function(rigid_ipc_download_gif_h)
   rigid_ipc_download_project(gif-h
    GIT_REPOSITORY https://github.com/charlietangora/gif-h.git
    GIT_TAG        606168840bff429a894dd470e36fb9956cbb07c0
  )
endfunction()

# Polysolve for Linear Solver Wrapper
function(rigid_ipc_download_polysolve)
  rigid_ipc_download_project(polysolve
    GIT_REPOSITORY https://github.com/polyfem/polysolve.git
    GIT_TAG        5d8f73476d79c79f3daf1bd831f1ad4ba2500711
  )
endfunction()

# IPC Toolkit
function(rigid_ipc_download_ipc_toolkit)
  rigid_ipc_download_project(ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG        614ee4584f53f422c17fd95995da1e6a9c06c41b
  )
endfunction()

# filib
function(rigid_ipc_download_filib)
  rigid_ipc_download_project(filib
    GIT_REPOSITORY https://github.com/txstc55/filib.git
    GIT_TAG        1b36c5047fd2ee27fcbb6d73deb25da7f77a6f70
  )
endfunction()

# SimpleBVH
function(rigid_ipc_download_simple_bvh)
  rigid_ipc_download_project(SimpleBVH
    GIT_REPOSITORY https://github.com/geometryprocessing/SimpleBVH.git
    GIT_TAG        15574502f6cb8039b0bfa4a85ccad04e09deaf05
  )
endfunction()

# tinygltf
function(rigid_ipc_download_tinygltf)
  rigid_ipc_download_project(tinygltf
    GIT_REPOSITORY https://github.com/syoyo/tinygltf.git
    GIT_TAG        v2.5.0
  )
endfunction()
