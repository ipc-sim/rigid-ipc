if(TARGET PolyFEM::polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'PolyFEM::polysolve'")

include(FetchContent)
FetchContent_Declare(
    polysolve
    GIT_REPOSITORY https://github.com/polyfem/polysolve.git
    GIT_TAG a94e9b8ed8302d4b479533c67419f31addb1e987
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(polysolve)

add_library(PolyFEM::polysolve ALIAS polysolve)
