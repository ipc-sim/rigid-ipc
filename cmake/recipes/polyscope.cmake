# Polyscope (https://github.com/nmwsharp/polyscope)
# License: MIT
if(TARGET polyscope::polyscope)
    return()
endif()

message(STATUS "Third-party: creating target 'polyscope::polyscope'")

include(CPM)
CPMAddPackage("gh:nmwsharp/polyscope@2.3.0")

add_library(polyscope::polyscope ALIAS polyscope)