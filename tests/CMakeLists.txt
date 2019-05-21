cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
    main.cpp

    # Test CCD
    collision_generator.cpp
    test_ccd.cpp
    test_impact_pruning.cpp
    test_collision_volume.cpp
    test_time_of_impact.cpp
    test_hash_grid.cpp

    # Test libarary integration
    integration/test_nlopt.cpp
    integration/test_osqp.cpp
    integration/test_mosek.cpp

    opt/test_ipopt_solver.cpp
    opt/test_ncp.cpp
    opt/test_ncp_displacements_opt.cpp
    opt/test_newton_solver.cpp
    opt/test_barriers.cpp
    opt/test_barrier_newton_solver.cpp

    opt/test_volume_constraint.cpp
    opt/test_barrier_constraint.cpp

    # opt/test_displacements_opt.cpp  # This file does not work yet
    opt/test_barrier_displacements_opt.cpp

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
