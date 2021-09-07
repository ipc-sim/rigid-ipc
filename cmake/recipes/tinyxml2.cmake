if(TARGET tinyxml2::tinyxml2)
    return()
endif()

message(STATUS "Third-party: creating target 'tinyxml2::tinyxml2'")

include(FetchContent)
FetchContent_Declare(
    tinyxml2
    GIT_REPOSITORY https://github.com/leethomason/tinyxml2.git
    GIT_TAG 9.0.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(tinyxml2)

add_library(tinyxml2::tinyxml2 ALIAS tinyxml2)
