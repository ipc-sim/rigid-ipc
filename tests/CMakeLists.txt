cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
    main.cpp
    collision_generator.cpp
    test_ccd.cpp
    test_impact_pruning.cpp
    test_collision_volume.cpp
    test_time_of_impact.cpp

    # opt/test_displacements_opt.cpp
    opt/test_newton_solver.cpp
    opt/test_ipopt_solver.cpp
    opt/test_osqp.cpp
    opt/test_mosek.cpp
)

add_executable(unit_tests ${test_sources})
target_include_directories(unit_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC ${PROJECT_NAME} catch warnings::all)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
