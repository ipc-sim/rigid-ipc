//#define _USE_MATH_DEFINES

//#include <cmath>
//#include <iostream>

//#include <catch2/catch.hpp>

// #include <igl/PI.h>
//#include <finitediff.hpp>

//#include <state.hpp>

//#include <Eigen/Geometry>

//#define NUM_ANGLES (5)

// using namespace ipc::rigid;
// using namespace opt;

// TEST_CASE("test the setup", "[opt][displacements][barrier]")
//{
//    State state;

//    // Load the test scene
//    state.load_scene(std::string(FIXTURES_DIR)
//        + "/tests/setup-falling-horizontal-edge.json");
//    REQUIRE(state.vertices.rows() == 4);
//    REQUIRE(state.displacements.rows() == 4);
//    REQUIRE(state.edges.rows() == 2);

//    state.constraint_function = GENERATE(ipc::rigid::ConstraintType::BARRIER);
//    state.getCollisionConstraint().update_collision_set = GENERATE(false,
//    true);

//    state.opt_method = ipc::rigid::OptimizationMethod::BARRIER_SOLVER;

//    state.reset_optimization_problem();

//    CHECK(state.particles_problem.num_vars == state.displacements.size());

//    CHECK(
//        state.particles_problem.x0.size() ==
//        state.particles_problem.num_vars);
//    CHECK(state.particles_problem.x_lower.size()
//        == state.particles_problem.num_vars);
//    CHECK(state.particles_problem.x_upper.size()
//        == state.particles_problem.num_vars);
//    CHECK(state.particles_problem.g_lower.size()
//        == state.particles_problem.num_constraints);
//    CHECK(state.particles_problem.g_upper.size()
//        == state.particles_problem.num_constraints);

//    state.particles_problem.x0.setRandom();

//    // Test ∇f
//    Eigen::VectorXd finite_grad(state.particles_problem.num_vars);
//    finite_gradient(state.particles_problem.x0,
//        state.particles_problem.func_f(), finite_grad);
//    Eigen::VectorXd analytic_grad
//        = state.particles_problem.eval_grad_f(state.particles_problem.x0);
//    CHECK(fd::compare_gradient(finite_grad, analytic_grad));

//    // Test ∇²f
//    Eigen::MatrixXd finite_hessian
//        = state.particles_problem.eval_hessian_f_approx(
//            state.particles_problem.x0);
//    Eigen::MatrixXd analytic_hessian
//        = state.particles_problem
//              .eval_hessian_f(state.particles_problem.x0)
//              .toDense();
//    CHECK(fd::compare_jacobian(finite_hessian, analytic_hessian));

//    // Test ∇g
//    Eigen::MatrixXd finite_jac
//        =
//        state.particles_problem.eval_jac_g_approx(state.particles_problem.x0);
//    Eigen::MatrixXd analytic_jac
//        = state.particles_problem.eval_jac_g(state.particles_problem.x0);
//    CHECK(fd::compare_jacobian(finite_jac, analytic_jac));

//    // Test ∇²g
//    Eigen::MatrixXd finite_hessian_gi(
//        state.particles_problem.num_vars, state.particles_problem.num_vars);
//    finite_hessian_gi.setOnes();
//    long i;

//    auto diff_i = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
//        auto jac = state.particles_problem.eval_jac_g(x);
//        return jac.row(i);
//    };

//    std::vector<Eigen::SparseMatrix<double>> analytic_hessian_g
//        = state.particles_problem.eval_hessian_g(state.particles_problem.x0);
//    for (i = 0; i < long(analytic_hessian_g.size()); i++) {
//        finite_jacobian(state.particles_problem.x0, diff_i,
//        finite_hessian_gi); CHECK(fd::compare_jacobian(
//            finite_hessian_gi, analytic_hessian_g[unsigned(i)]));
//    }
//}

// TEST_CASE("two rotating edges", "[opt][displacements][barrier]")
//{
//    State state;
//    state.load_scene(std::string(FIXTURES_DIR)
//        + "/two-edges/falling-edge-00°-same-length.json");

//    REQUIRE(state.vertices.rows() == 4);
//    REQUIRE(state.displacements.rows() == 4);
//    REQUIRE(state.edges.rows() == 2);

//    state.constraint_function = ipc::rigid::ConstraintType::BARRIER;

//    state.getCollisionConstraint().update_collision_set = GENERATE(false,
//    true); state.opt_method = ipc::rigid::OptimizationMethod::BARRIER_SOLVER;

//    double theta1 = 2 * igl::PI / NUM_ANGLES * GENERATE(range(0, NUM_ANGLES));
//    double theta2 = 2 * igl::PI / NUM_ANGLES * GENERATE(range(0, NUM_ANGLES));

//    for (int i = 0; i < state.edges.rows(); i++) {
//        Eigen::Rotation2D<double> R(i == 0 ? theta1 : theta2);
//        Eigen::RowVector2d center = (state.vertices.row(state.edges(i, 1))
//                                        + state.vertices.row(state.edges(i,
//                                        0)))
//            / 2;
//        state.vertices.row(state.edges(i, 0))
//            = (R * (state.vertices.row(state.edges(i, 0)) -
//            center).transpose())
//                  .transpose()
//            + center;
//        state.vertices.row(state.edges(i, 1))
//            = (R * (state.vertices.row(state.edges(i, 1)) -
//            center).transpose())
//                  .transpose()
//            + center;
//    }

//    state.barrier_solver.min_barrier_epsilon = 1e-3;
//    state.barrier_solver.inner_solver_type = BarrierInnerSolver::NEWTON;
//    dynamic_cast<NewtonSolver&>(state.barrier_solver.get_inner_solver())
//        .min_step_length
//        = 1e-4;
//    dynamic_cast<NewtonSolver&>(state.barrier_solver.get_inner_solver())
//        .absolute_tolerance
//        = 1e-3;

//    state.optimize_displacements("");
//    CHECK(state.opt_results.success);
//}

// TEST_CASE("corner case", "[opt][displacements][barrier]")
//{
//    State state;
//    state.load_scene(std::string(FIXTURES_DIR) + "/corner-case.json");

//    REQUIRE(state.vertices.rows() == 5);
//    REQUIRE(state.displacements.rows() == 5);
//    REQUIRE(state.edges.rows() == 3);

//    state.constraint_function = ipc::rigid::ConstraintType::BARRIER;

//    state.getCollisionConstraint().update_collision_set = GENERATE(false,
//    true); state.opt_method = ipc::rigid::OptimizationMethod::BARRIER_SOLVER;

//    state.barrier_solver.min_barrier_epsilon = 1e-3;
//    state.barrier_solver.inner_solver_type = BarrierInnerSolver::NEWTON;
//    dynamic_cast<NewtonSolver&>(state.barrier_solver.get_inner_solver())
//        .min_step_length
//        = 1e-4;
//    dynamic_cast<NewtonSolver&>(state.barrier_solver.get_inner_solver())
//        .absolute_tolerance
//        = 1e-3;

//    state.optimize_displacements("");
//    CHECK(state.opt_results.success);
//}
