################################################################################
# Tests
################################################################################

add_executable(rigid_ipc_tests
  main.cpp

  barrier/test_barriers.cpp

  # Test CCD
  ccd/collision_generator.cpp
  ccd/test_edge_vertex_ccd.cpp
  ccd/test_time_of_impact.cpp
  # ccd/test_hash_grid.cpp

  interval/test_interval.cpp
  interval/test_interval_root_finder.cpp
  ccd/test_rigid_body_time_of_impact.cpp
  ccd/test_rigid_body_hash_grid.cpp

  solvers/test_newton_solver.cpp
  solvers/test_barrier_newton_solver.cpp
  solvers/test_barrier_displacements_opt.cpp

  opt/test_distance_barrier_constraint.cpp

  physics/test_mass.cpp
  physics/test_pose.cpp
  physics/test_rigid_body.cpp
  physics/test_rigid_body_system.cpp
  physics/test_rigid_body_problem.cpp

  io/test_serialize_json.cpp
  io/test_read_rb_scene.cpp

  geometry/test_distance.cpp
  geometry/test_intersection.cpp

  utils/test_sinc.cpp
)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(rigid_ipc_tests PUBLIC ipc::rigid)

include(rigid_ipc_warnings)
target_link_libraries(rigid_ipc_tests PRIVATE ipc::rigid::warnings)

include(catch2)
target_link_libraries(rigid_ipc_tests PUBLIC Catch2::Catch2)

include(finite_diff)
target_link_libraries(rigid_ipc_tests PUBLIC finitediff::finitediff)

################################################################################
# Compiler options
################################################################################

target_compile_definitions(rigid_ipc_tests PUBLIC CATCH_CONFIG_ENABLE_BENCHMARKING)

################################################################################
# Register tests
################################################################################

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(rigid_ipc_tests)
