#include <catch2/catch.hpp>

#include <Eigen/Eigenvalues>
#include <utils/eigen_ext.hpp>

#include <opt/ad_hoc_problem.hpp>
#include <solvers/newton_solver.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Simple tests of Newton's Method", "[opt][newtons_method]")
{
    int num_vars = GENERATE(1, 10, 100), num_constraints = 0;
    // Setup solver
    NewtonSolver solver;
    solver.free_dof = Eigen::VectorXi::LinSpaced(num_vars, 0, num_vars - 1);
    // Setup problem
    // -----------------------------------------------------------------
    AdHocProblem problem(num_vars, num_constraints);
    problem.x0.setRandom();

    problem.f = [](const Eigen::VectorXd& x) { return x.squaredNorm() / 2.0; };
    problem.grad_f = [](const Eigen::VectorXd& x) { return x; };
    problem.hessian_f = [](const Eigen::VectorXd& x) {
        return Eigen::MatrixXd::Identity(x.rows(), x.rows());
    };

    REQUIRE(problem.validate_problem());

    OptimizationResults results = solver.solve(problem);
    REQUIRE(results.success);
    CHECK(results.x.squaredNorm() == Approx(0).margin(1e-6));
    CHECK(results.minf == Approx(0).margin(1e-6));
}

TEST_CASE("Test Newton direction solve", "[opt][newtons_method][newton_dir]")
{
    int num_vars = 1000;
    Eigen::VectorXd x(num_vars);
    x.setRandom();
    // f = x^2
    Eigen::VectorXd gradient = 2 * x;
    Eigen::SparseMatrix<double> hessian
        = Eigen::SparseDiagonal<double>(2 * Eigen::VectorXd::Ones(num_vars));
    Eigen::VectorXd delta_x;
    ccd::opt::NewtonSolver solver;
    solver.compute_direction(gradient, hessian, delta_x);
    CHECK((x + delta_x).squaredNorm() == Approx(0.0));
}

TEST_CASE("Test making a matrix SPD", "[opt][make_spd]")
{
    Eigen::SparseMatrix<double> A
        = Eigen::MatrixXd::Random(100, 100).sparseView();
    double mu = ccd::opt::make_matrix_positive_definite(A);
    CAPTURE(mu);
    auto eig_vals = Eigen::MatrixXd(A).eigenvalues();
    for (int i = 0; i < eig_vals.size(); i++) {
        CHECK(eig_vals(i).real() >= Approx(0.0).margin(1e-12));
    }
}
