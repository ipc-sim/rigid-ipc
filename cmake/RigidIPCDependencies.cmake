# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(RigidIPCDownloadExternal)

################################################################################
# Required libraries
################################################################################

# format strings
if(NOT TARGET fmt::fmt)
  rigid_ipc_download_fmt()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/fmt)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
  rigid_ipc_download_spdlog()
  add_library(spdlog INTERFACE)
  add_library(spdlog::spdlog ALIAS spdlog)
  target_include_directories(spdlog SYSTEM INTERFACE ${RIGID_IPC_EXTERNAL}/spdlog/include)
  target_compile_definitions(spdlog INTERFACE -DSPDLOG_FMT_EXTERNAL)
  target_link_libraries(spdlog INTERFACE fmt::fmt)
endif()

# libigl
if(NOT TARGET igl::core)
  rigid_ipc_download_libigl()
  # Import libigl targets
  list(APPEND CMAKE_MODULE_PATH "${RIGID_IPC_EXTERNAL}/libigl/cmake")
  include(libigl)

  if(NOT RIGID_IPC_WITH_OPENGL AND NOT TARGET stb_image)
    # Download this myself for the software renderer
    igl_download_stb()
    add_subdirectory(${LIBIGL_EXTERNAL}/stb stb_image)
  endif()
endif()

# json
if(NOT TARGET nlohmann::json)
  rigid_ipc_download_json()
  option(JSON_BuildTests "" OFF)
  option(JSON_MultipleHeaders "" ON)
  add_subdirectory(${RIGID_IPC_EXTERNAL}/json json)
  add_library(nlohmann::json ALIAS nlohmann_json)
endif()

if(NOT TARGET CLI11::CLI11)
  rigid_ipc_download_cli11()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/cli11)
endif()

# finite-diff
if(NOT TARGET FiniteDiff::FiniteDiff)
  rigid_ipc_download_finite_diff()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/finite-diff EXCLUDE_FROM_ALL)
  add_library(FiniteDiff::FiniteDiff ALIAS FiniteDiff)
endif()

# TBB
if(NOT TARGET TBB::tbb)
  rigid_ipc_download_tbb()
  set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
  set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
  set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${RIGID_IPC_EXTERNAL}/tbb EXCLUDE_FROM_ALL)
  add_library(TBB::tbb ALIAS tbb_static)
endif()

# Tight Inclusion for MSCCD
if(NOT TARGET TightInclusion)
  rigid_ipc_download_tight_inclusion()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/Tight-Inclusion)
  add_library(TightInclusion ALIAS tight_inclusion)
endif()

# XML for Bullet/MuJoCo comparisons
if(NOT TARGET tinyxml2::tinyxml2)
  rigid_ipc_download_tinyxml2()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/tinyxml2)
  add_library(tinyxml2::tinyxml2 ALIAS tinyxml2)
endif()

# GIF writer for the viewer
if(NOT TARGET gif_h::gif_h)
  rigid_ipc_download_gif_h()
  add_library(gif_h INTERFACE)
  target_include_directories(gif_h INTERFACE ${RIGID_IPC_EXTERNAL}/gif-h)
  add_library(gif_h::gif_h ALIAS gif_h)
endif()

# Polysolve for Linear Solver Wrapper
if(NOT TARGET PolyFEM::polysolve)
  rigid_ipc_download_polysolve()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/polysolve)
  add_library(PolyFEM::polysolve ALIAS polysolve)
endif()

# IPC Toolkit
if(NOT TARGET ipc::toolkit)
  rigid_ipc_download_ipc_toolkit()
  set(IPC_TOOLKIT_BUILD_UNIT_TESTS OFF CACHE BOOL " " FORCE)
  add_subdirectory(${RIGID_IPC_EXTERNAL}/ipc-toolkit)
  add_library(ipc::toolkit ALIAS IPCToolkit)
endif()

# filib
if(NOT TARGET filib::filib)
  rigid_ipc_download_filib()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/filib)
  add_library(filib::filib ALIAS filib)
endif()

# SimpleBVH
if(NOT TARGET SimpleBVH::BVH_lib)
  rigid_ipc_download_simple_bvh()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/SimpleBVH)
  add_library(SimpleBVH::BVH_lib ALIAS BVH_lib)
endif()

# tinygltf
if(NOT TARGET tinygltf::tinygltf)
  rigid_ipc_download_tinygltf()
  add_library(tinygltf INTERFACE)
  target_include_directories(tinygltf INTERFACE ${RIGID_IPC_EXTERNAL}/tinygltf)
  add_library(tinygltf::tinygltf ALIAS tinygltf)
endif()

# GHC Filesystem
if(NOT TARGET ghc::filesystem)
  rigid_ipc_download_filesystem()
  add_subdirectory(${RIGID_IPC_EXTERNAL}/filesystem)
  add_library(ghc::filesystem ALIAS ghc_filesystem)
endif()
