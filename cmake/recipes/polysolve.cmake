if(TARGET PolyFEM::polysolve)
    return()
endif()

message(STATUS "Third-party: creating target 'PolyFEM::polysolve'")

include(FetchContent)
FetchContent_Declare(
    polysolve
    GIT_REPOSITORY https://github.com/polyfem/polysolve.git
    GIT_TAG 72e5eaca17b1ae975fa5a7149627a17e6b13cf80
    GIT_SHALLOW FALSE
)
FetchContent_MakeAvailable(polysolve)

add_library(PolyFEM::polysolve ALIAS polysolve)
