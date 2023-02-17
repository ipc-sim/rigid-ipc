if(TARGET software_renderer::software_renderer)
    return()
endif()

message(STATUS "Third-party: creating target 'software_renderer::software_renderer'")

include(FetchContent)
FetchContent_Declare(
    software_renderer
    GIT_REPOSITORY https://github.com/zfergus/software-renderer.git
    GIT_TAG 2453b2b05b42e95ff38d1b6f91c9857390087a69
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(software_renderer)
