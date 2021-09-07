if(TARGET tinygltf::tinygltf)
    return()
endif()

message(STATUS "Third-party: creating target 'tinygltf::tinygltf'")

include(FetchContent)
FetchContent_Declare(
    tinygltf
    GIT_REPOSITORY https://github.com/syoyo/tinygltf.git
    GIT_TAG v2.5.0
    GIT_SHALLOW TRUE
)

FetchContent_GetProperties(tinygltf)
if(NOT tinygltf_POPULATED)
    FetchContent_Populate(tinygltf)
endif()

add_library(tinygltf INTERFACE)
target_include_directories(tinygltf INTERFACE ${tinygltf_SOURCE_DIR})
add_library(tinygltf::tinygltf ALIAS tinygltf)
