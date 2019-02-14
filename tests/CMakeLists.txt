cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
    main.cpp
    test_ccd.cpp
    test_impact_pruning.cpp
    test_collision_volume.cpp
    test_autogen.cpp
    test_compute_edge_vertex_time_of_impact.cpp
)

add_executable(unit_tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC ${PROJECT_NAME} catch warnings::all)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)



if(BUILD_COVERAGE)
    include(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    SETUP_TARGET_FOR_COVERAGE_LCOV(
        NAME ${PROJECT_NAME}_coverage EXECUTABLE unit_tests DEPENDS ${PROJECT_NAME})
# might work too ... but just as bad
#    SETUP_TARGET_FOR_COVERAGE_LCOV(
#        NAME ${PROJECT_NAME}_coverage EXECUTABLE ctest)
endif()
