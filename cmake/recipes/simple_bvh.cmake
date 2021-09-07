if(TARGET SimpleBVH::BVH_lib)
    return()
endif()

message(STATUS "Third-party: creating target 'SimpleBVH::BVH_lib'")

include(FetchContent)
FetchContent_Declare(
    SimpleBVH
    GIT_REPOSITORY https://github.com/geometryprocessing/SimpleBVH.git
    GIT_TAG 15574502f6cb8039b0bfa4a85ccad04e09deaf05
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(SimpleBVH)

add_library(SimpleBVH::BVH_lib ALIAS BVH_lib)
