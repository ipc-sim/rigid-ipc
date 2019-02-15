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
