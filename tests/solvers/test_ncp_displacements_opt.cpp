//#include <catch2/catch.hpp>

//#include <iostream>
//#include <state.hpp>

////-----------------
//// Tests
//// ---------------------------------------------------

// const int NUM_EDGES = 2;
// const int NUM_VERTICES = 4;

// TEST_CASE("NCP Displacements", "[opt][NCP][NCP-Displacements]")
//{
//    ccd::State state;

//    state.opt_method = ccd::OptimizationMethod::NCP;
//    state.ncp_solver.max_iterations = 200;
//    state.ncp_solver.convergence_tolerance = 1E-8;
//    state.ncp_solver.check_convergence = false;
//    state.ncp_solver.solve_for_active_cstr = false;
//    state.particles_problem.use_mass_matrix = false;

//    state.constraint_function = ccd::ConstraintType::VOLUME;
//    state.volume_constraint.volume_epsilon = 1E-8;

//    Eigen::MatrixX2d vertices(NUM_VERTICES, 2);
//    Eigen::MatrixX2i edges(NUM_EDGES, 2);
//    Eigen::MatrixX2d displacements(NUM_VERTICES, 2);

//    edges.row(0) << 0, 1;
//    edges.row(1) << 2, 3;

//    // horizontal fixed edge of length 1.0
//    vertices.row(0) << -0.5, 0.0;
//    vertices.row(1) << 0.5, 0.0;

//    displacements.row(0) << 0.0, 0.0;
//    displacements.row(1) << 0.0, 0.0;

//    Eigen::MatrixX2d expected(NUM_VERTICES, 2);

//    SECTION("Vertical Displ")
//    {
//        vertices.row(2) << 0.0, 0.5;
//        vertices.row(3) << 0.0, 1.0;
//        displacements.row(2) << 0.0, -0.6;
//        displacements.row(3) << 0.0, -0.6;

//        expected << 2.42861e-17, -0.0333333, -2.42861e-17, -0.0333333, 0,
//            -0.533333, 0, -0.6;
//    }

//    SECTION("Vertical Displ Long")
//    {
//        vertices.row(2) << 0.0, 0.5;
//        vertices.row(3) << 0.0, 1.0;
//        displacements.row(2) << 0.0, -1.0;
//        displacements.row(3) << 0.0, -1.0;

//        expected << 0.000626412, 1.23514e-08, -0.000626412, -1.23514e-08,
//            -5.01155e-07, -0.890427, -6.67193e-08, -1.00049;
//    }

//    state.load_scene(vertices, edges, displacements);
//    state.optimize_displacements();

//    Eigen::MatrixX2d u_ = state.opt_results.x;
//    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols,
//        ", ", ", ", "", "", " << ", ";");

//    CHECK((u_ - expected).squaredNorm() < 1E-4);
//}
