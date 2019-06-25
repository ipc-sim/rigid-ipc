cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(unit_tests
    main.cpp

    # Test CCD
    ccd/collision_generator.cpp
    ccd/test_ccd.cpp
    ccd/test_impact_pruning.cpp
    ccd/test_time_of_impact.cpp
    ccd/test_hash_grid.cpp

    # Test libarary integration
    integration/test_nlopt.cpp
    integration/test_osqp.cpp
    integration/test_mosek.cpp
    integration/test_lbfgspp.cpp

    solvers/test_ipopt_solver.cpp
    solvers/test_ncp.cpp
    solvers/test_ncp_displacements_opt.cpp
    solvers/test_newton_solver.cpp
    solvers/test_barrier_newton_solver.cpp
    solvers/test_barrier_displacements_opt.cpp
    solvers/test_line_search.cpp

    opt/test_barriers.cpp
    opt/test_volume_constraint.cpp
    opt/test_barrier_constraint.cpp
    opt/test_rigid_body_problem.cpp

    physics/test_rigid_body.cpp
    physics/test_rigid_body_system.cpp

)

target_link_libraries(unit_tests PUBLIC ${PROJECT_NAME} Catch2::Catch2 warnings::all)

if(CODE_COVERAGE)
    target_code_coverage(unit_tests AUTO ALL)
endif()

foreach(source IN ITEMS ${test_sources})
   source_group("tests" FILES "${source}")
endforeach()

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
