if(TARGET tight_inclusion::tight_inclusion)
    return()
endif()

message(STATUS "Third-party: creating target 'tight_inclusion::tight_inclusion'")

include(FetchContent)
FetchContent_Declare(
    tight_inclusion
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Tight-Inclusion.git
    GIT_TAG b913416aecf1bbb8a51726039dcb94f395849cab
    GIT_SHALLOW FALSE
)

FetchContent_MakeAvailable(tight_inclusion)

add_library(tight_inclusion::tight_inclusion ALIAS tight_inclusion)
