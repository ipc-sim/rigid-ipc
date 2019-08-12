#include <catch2/catch.hpp>

#include <Eigen/Eigenvalues>
#include <utils/eigen_ext.hpp>

#include <opt/ad_hoc_problem.hpp>
#include <solvers/newton_solver.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("Simple tests of Newton's Method", "[opt][newtons_method]")
{
    int num_vars = GENERATE(1, 10, 100);

    NewtonSolver solver;
    solver.init_free_dof(Eigen::VectorXb::Zero(num_vars));

    // Setup problem
    // -----------------------------------------------------------------
    class AdHocProblem : public virtual IBarrierProblem {
    public:
        AdHocProblem(int num_vars)
        {
            num_vars_ = num_vars;
            x0.resize(num_vars);
            x0.setRandom();
            is_dof_fixed_ = Eigen::VectorXb::Zero(num_vars);
        }
        double eval_f(const Eigen::VectorXd& x) override
        {
            return x.squaredNorm() / 2.0;
        }
        double eval_f_(const Eigen::VectorXd& x,
            const bool /*update_constraint_set*/) override
        {
            return eval_f(x);
        }

        Eigen::VectorXd eval_grad_f(const Eigen::VectorXd& x) override
        {
            return x;
        }
        Eigen::SparseMatrix<double> eval_hessian_f(
            const Eigen::VectorXd& x) override
        {
            return Eigen::MatrixXd::Identity(x.rows(), x.rows()).sparseView();
        }
        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian) override
        {
            f_uk = eval_f(x);
            f_uk_jacobian = eval_grad_f(x);
        }
        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override
        {
            f_uk = eval_f(x);
            f_uk_jacobian = eval_grad_f(x);
            f_uk_hessian = eval_hessian_f(x);
        }
        double get_barrier_epsilon() override { return 1.0; }
        const Eigen::VectorXd& starting_point() override { return x0; }
        const int& num_vars() override { return num_vars_; }
        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }

        int num_vars_;
        Eigen::VectorXb is_dof_fixed_;
        Eigen::VectorXd x0;
    };

    AdHocProblem problem(num_vars);

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
