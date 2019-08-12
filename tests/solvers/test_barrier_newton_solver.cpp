#include <catch2/catch.hpp>
#include <cmath>

#include <autodiff/finitediff.hpp>
#include <opt/ad_hoc_problem.hpp>
#include <solvers/barrier_solver.hpp>
#include <solvers/newton_solver.hpp>

#include <logger.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE(
    "Check barrier problem derivatives", "[opt][barrier][barrier_problem]")
{
    // TODO: Generate a random problem and test the derivatives
    int num_vars = GENERATE(1, 2, 5, 10), num_constraints = 1;
    class AdHocProblem : public virtual IBarrierGeneralProblem {
    public:
        AdHocProblem(int num_vars, int num_constraints, double epsilon)
            : num_vars_(num_vars)
            , num_constraints_(num_constraints)
            , eps(epsilon)
        {
        }
        double eval_f(const Eigen::VectorXd& x) override
        {
            return x.squaredNorm() / 2;
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
            Eigen::VectorXd& f_uk_jacobian,
            Eigen::SparseMatrix<double>& f_uk_hessian) override
        {
            f_uk = eval_f(x);
            f_uk_jacobian = eval_grad_f(x);
            f_uk_hessian = eval_hessian_f(x);
        }
        void eval_f_and_fdiff(const Eigen::VectorXd& x,
            double& f_uk,
            Eigen::VectorXd& f_uk_jacobian) override
        {
            f_uk = eval_f(x);
            f_uk_jacobian = eval_grad_f(x);
        }
        Eigen::VectorXd eval_g(const Eigen::VectorXd& x) override
        {
            Eigen::VectorXd gx(num_constraints_);
            gx.setConstant(eval_f(x));
            return gx;
        }
        Eigen::VectorXd eval_g_(const Eigen::VectorXd& x,
            const bool /*update_constraint_set*/) override
        {
            return eval_g(x);
        }
        Eigen::MatrixXd eval_jac_g(const Eigen::VectorXd& x) override
        {
            return eval_grad_f(x).transpose();
        }
        std::vector<Eigen::SparseMatrix<double>> eval_hessian_g(
            const Eigen::VectorXd& x) override
        {
            std::vector<Eigen::SparseMatrix<double>> hessian;
            for (long i = 0; i < num_constraints_; i++) {
                hessian.push_back(eval_hessian_f(x));
            }
            return hessian;
        }
        void eval_g_and_gdiff(const Eigen::VectorXd& x,
            Eigen::VectorXd& gx,
            Eigen::MatrixXd& gx_jacobian,
            std::vector<Eigen::SparseMatrix<double>>& gx_hessian) override
        {
            gx = eval_g(x);
            gx_jacobian = eval_jac_g(x);
            gx_hessian = eval_hessian_g(x);
        }

        double get_barrier_epsilon() override { return eps; }
        void set_barrier_epsilon(const double epsilon) override
        {
            eps = epsilon;
        }
        const int& num_vars() override { return num_vars_; }
        const int& num_constraints() override { return num_constraints_; }
        const Eigen::VectorXb& is_dof_fixed() override { return is_dof_fixed_; }
        const Eigen::VectorXd& starting_point() override { return x0; }

        int num_vars_;
        int num_constraints_;
        double eps;
        Eigen::VectorXd x0;
        Eigen::VectorXb is_dof_fixed_;
    };

    double epsilon = GENERATE(1.0, 0.5, 1e-1, 5e-2);
    AdHocProblem problem(num_vars, num_constraints, epsilon);

    BarrierProblem barrier_problem(problem);

    Eigen::VectorXd x(num_vars);
    x.setConstant(GENERATE(-1.0, 1.0) * GENERATE(1e-1, 1.0, 2.0 - 1e-3, 4.0));

    // If the function evaluates to infinity then the finite differences will
    // not work. I assume in the definition of the barrier gradient that
    // d/dx ∞ = 0.
    if (!std::isinf(barrier_problem.eval_f(x))) {
        // Use a higher order finite difference method because the function near
        // the boundary becomes very non-linear. This problem worsens as the ϵ
        // of the boundary gets smaller.

        // Test ∇f
        Eigen::VectorXd finite_grad(barrier_problem.num_vars());
        finite_grad = ccd::opt::eval_grad_f_approx(barrier_problem, x);
        Eigen::VectorXd analytic_grad = barrier_problem.eval_grad_f(x);
        CHECK(compare_gradient(finite_grad, analytic_grad));

        // Test ∇²f
        Eigen::MatrixXd finite_hessian = eval_hess_f_approx(barrier_problem, x);

        Eigen::MatrixXd analytic_hessian
            = barrier_problem.eval_hessian_f(x).toDense();
        CHECK(compare_jacobian(finite_hessian, analytic_hessian));

        CAPTURE(x, problem.eval_g(x), epsilon, barrier_problem.eval_f(x),
            finite_grad, analytic_grad, finite_hessian, analytic_hessian);
    }
}
