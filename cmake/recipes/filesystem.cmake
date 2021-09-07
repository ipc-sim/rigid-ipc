if(TARGET ghc::filesystem)
    return()
endif()

message(STATUS "Third-party: creating target 'ghc::filesystem'")

include(FetchContent)
FetchContent_Declare(
    filesystem
    GIT_REPOSITORY https://github.com/gulrak/filesystem.git
    GIT_TAG v1.5.8
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(filesystem)
add_library(ghc::filesystem ALIAS ghc_filesystem)
