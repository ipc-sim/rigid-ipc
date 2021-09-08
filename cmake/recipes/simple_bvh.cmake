if(TARGET simple_bvh::simple_bvh)
    return()
endif()

message(STATUS "Third-party: creating target 'simple_bvh::simple_bvh'")

include(FetchContent)
FetchContent_Declare(
    SimpleBVH
    GIT_REPOSITORY https://github.com/ipc-sim/SimpleBVH.git
    GIT_TAG 6d4b162cc87f156cf9a538a447ef17b94233d60c
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(SimpleBVH)
