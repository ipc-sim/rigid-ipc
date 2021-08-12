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
        GIT_TAG        8.0.1
    )
endfunction()

function(rigid_ipc_download_spdlog)
     rigid_ipc_download_project(spdlog
       GIT_REPOSITORY https://github.com/gabime/spdlog.git
       GIT_TAG        v1.9.1
    )
endfunction()

function(rigid_ipc_download_libigl)
     rigid_ipc_download_project(libigl
       GIT_REPOSITORY https://github.com/libigl/libigl.git
       GIT_TAG        v2.3.0
    )
endfunction()

function(rigid_ipc_download_json)
     rigid_ipc_download_project(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.9.1
    )
endfunction()

function(rigid_ipc_download_cli11)
     rigid_ipc_download_project(cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG        v2.0.0
    )
endfunction()

function(rigid_ipc_download_finite_diff)
     rigid_ipc_download_project(finite-diff
        GIT_REPOSITORY https://github.com/zfergus/finite-diff.git
        GIT_TAG        f35375d2db00618d19ed7e4d2ce505288006f403
    )
endfunction()

function(rigid_ipc_download_tbb)
   rigid_ipc_download_project(tbb
    GIT_REPOSITORY https://github.com/wjakob/tbb.git
    GIT_TAG        9e219e24fe223b299783200f217e9d27790a87b0
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
    GIT_TAG        9.0.0
  )
endfunction()

# GIF writer
function(rigid_ipc_download_gif_h)
   rigid_ipc_download_project(gif-h
    GIT_REPOSITORY https://github.com/charlietangora/gif-h.git
    GIT_TAG        3d2657b9ad40aac9fd6f75ad079335856e94d664
  )
endfunction()

# Polysolve for Linear Solver Wrapper
function(rigid_ipc_download_polysolve)
  rigid_ipc_download_project(polysolve
    GIT_REPOSITORY https://github.com/polyfem/polysolve.git
    GIT_TAG        a94e9b8ed8302d4b479533c67419f31addb1e987
  )
endfunction()

# IPC Toolkit
function(rigid_ipc_download_ipc_toolkit)
  rigid_ipc_download_project(ipc-toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG        24056ccb2ca0a03bdef8141bc5011c41547f06b5
  )
endfunction()

# filib
function(rigid_ipc_download_filib)
  rigid_ipc_download_project(filib
    GIT_REPOSITORY https://github.com/txstc55/filib.git
    GIT_TAG        31c801fc45b545a7809911b193498efdb4bd930d
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

# Filesystem library for C++11 and C++14
function(rigid_ipc_download_filesystem)
  rigid_ipc_download_project(filesystem
    GIT_REPOSITORY https://github.com/gulrak/filesystem.git
    GIT_TAG        v1.5.8
  )
endfunction()
