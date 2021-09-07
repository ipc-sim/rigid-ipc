if(TARGET software_renderer::software_renderer)
    return()
endif()

message(STATUS "Third-party: creating target 'software_renderer::software_renderer'")

include(FetchContent)
FetchContent_Declare(
    software_renderer
    GIT_REPOSITORY https://github.com/zfergus/software-renderer.git
    GIT_TAG cdfb409d02a1b3b55e916c549aedbff7e5bdb5e7
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(software_renderer)
