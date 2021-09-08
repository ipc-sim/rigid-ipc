if(TARGET gif_h::gif_h)
    return()
endif()

message(STATUS "Third-party: creating target 'gif_h::gif_h'")

include(FetchContent)
FetchContent_Declare(
    gif_h
    GIT_REPOSITORY https://github.com/charlietangora/gif-h.git
    GIT_TAG 3d2657b9ad40aac9fd6f75ad079335856e94d664
    GIT_SHALLOW FALSE
)

FetchContent_GetProperties(gif_h)
if(NOT gif_h_POPULATED)
    FetchContent_Populate(gif_h)
endif()

add_library(gif_h INTERFACE)
target_include_directories(gif_h INTERFACE ${gif_h_SOURCE_DIR})
add_library(gif_h::gif_h ALIAS gif_h)
