// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include "barrier_newton_solver.hpp"

#include <opt/barrier.hpp>
#include <solvers/newtons_method.hpp>
#include <utils/eigen_ext.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    BarrierProblem::BarrierProblem(
        OptimizationProblem& general_problem, const double epsilon)
        : general_problem(&general_problem)
        , epsilon(epsilon)
    {
        num_vars = general_problem.num_vars;
        num_constraints = general_problem.num_constraints;
        x0 = general_problem.x0;
        fixed_dof = general_problem.fixed_dof;
    }

    BarrierProblem::~BarrierProblem() {}

    Eigen::VectorXd BarrierProblem::barrier(const Eigen::VectorXd x)
    {
        const double eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier(xi, eps); });
    }

    Eigen::VectorXd BarrierProblem::barrier_gradient(const Eigen::VectorXd x)
    {
        const double eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier_gradient(xi, eps); });
    }

    Eigen::VectorXd BarrierProblem::barrier_hessian(const Eigen::VectorXd x)
    {
        const double eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier_hessian(xi, eps); });
    }

    void BarrierProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk, Eigen::VectorXd& f_uk_gradient,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        general_problem->eval_f_and_fdiff(x, f_uk, f_uk_gradient, f_uk_hessian);

        Eigen::VectorXd gx;
        Eigen::MatrixXd dgx;
        std::vector<Eigen::SparseMatrix<double>> ddgx;
        general_problem->eval_g_and_gdiff(x, gx, dgx, ddgx);

        // Check the sizes of the x's match
        assert(x.size() == general_problem->x_lower.size());
        assert(x.size() == general_problem->x_upper.size());

        // Sum the barrier constraints into the objective
        f_uk += gx.sum() + barrier(x - general_problem->x_lower).sum()
            + barrier(-x + general_problem->x_upper).sum();

        // Sum the barrier constraint gradients into the objective gradient
        f_uk_gradient += dgx.colwise().sum().transpose()
            + barrier_gradient(x - general_problem->x_lower)
            - barrier_gradient(-x + general_problem->x_upper);

        // Sum the barrier constraint hessian into the objective hessian
        PROFILE(
            // clang-format off
            for (const auto& ddgx_i : ddgx) {
                f_uk_hessian += ddgx_i;
            },
            // clang-format on
            ProfiledPoint::SUMMING_HESSIAN)
    }

    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {
        // Check the sizes of the x's match
        assert(x.size() == general_problem->x_lower.size());
        assert(x.size() == general_problem->x_upper.size());

        double val = general_problem->eval_f(x)
            + barrier(x - general_problem->x_lower).sum()
            + barrier(-x + general_problem->x_upper).sum();

        Eigen::VectorXd gx = general_problem->eval_g(x);
        val += gx.sum();

        return val;
    }

    Eigen::VectorXd BarrierProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd grad = general_problem->eval_grad_f(x);

        // Add constraint functions barrier(g(x))
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);
        grad += dgx.colwise().sum().transpose();

        // Add value constraint barriers
        // ∇ ∑ ϕ(x_i) = ∑ (∇ ϕ(x_i)) = ∑ [0 ... ϕ'(x_i) ... 0]^T
        //            = [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
        // Check the sizes of the x's match
        assert(x.size() == general_problem->x_lower.size());
        assert(x.size() == general_problem->x_upper.size());
        grad += barrier_gradient(x - general_problem->x_lower)
            - barrier_gradient(-x + general_problem->x_upper);
        return grad;
    }

    Eigen::MatrixXd BarrierProblem::eval_hessian_f(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd hessian = general_problem->eval_hessian_f(x);

        // Eigen::VectorXd gx = general_problem->eval_g(x);
        // Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);
        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x);

        PROFILE(
            // clang-format off
            Eigen::SparseMatrix<double> sum_ddgx(
                hessian.rows(), hessian.cols());
            for (const auto& ddgx_i : ddgx) {
                sum_ddgx += ddgx_i;
            }
            hessian += sum_ddgx;
            // clang-format on
            , ProfiledPoint::SUMMING_HESSIAN)

        // ∇ [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
        // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
        // Check the sizes of the x's match
        assert(x.size() == general_problem->x_lower.size());
        assert(x.size() == general_problem->x_upper.size());
        hessian.diagonal() += barrier_hessian(x - general_problem->x_lower)
            + barrier_hessian(-x + general_problem->x_upper);
        return hessian;
    }

    Eigen::SparseMatrix<double> BarrierProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        Eigen::SparseMatrix<double> hessian
            = general_problem->eval_hessian_f_sparse(x); // mass matrix

        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x); // hessian of constraints

        PROFILE(
            // clang-format off
            for (const auto& ddgx_i : ddgx) {
                hessian += ddgx_i;
            },
            // clang-format on
            ProfiledPoint::SUMMING_HESSIAN)

        // ∇ [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
        // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
        assert(x.size() == general_problem->x_lower.size());
        assert(x.size() == general_problem->x_upper.size());
        hessian += Eigen::SparseDiagonal<double>(
            barrier_hessian(x - general_problem->x_lower)
            + barrier_hessian(-x + general_problem->x_upper));
        return hessian;
    }

    void BarrierProblem::enable_line_search_mode(const Eigen::VectorXd& max_x)
    {
        general_problem->enable_line_search_mode(max_x);
    }

    void BarrierProblem::disable_line_search_mode()
    {
        general_problem->disable_line_search_mode();
    }

    bool BarrierProblem::eval_intermediate_callback(const Eigen::VectorXd& x)
    {
        return general_problem->eval_intermediate_callback(x);
    }

    BarrierNewtonSolver::BarrierNewtonSolver()
        : barrier_constraint(nullptr)
        , min_barrier_epsilon(1e-5)
        , absolute_tolerance(1e-5)
        , line_search_tolerance(1e-12)
        , max_iterations(3000)
    {
    }

    BarrierNewtonSolver::~BarrierNewtonSolver() {}

    OptimizationResults BarrierNewtonSolver::solve(
        OptimizationProblem& general_problem)
    {
        assert(barrier_constraint != nullptr);

        BarrierProblem barrier_problem(
            general_problem, barrier_constraint->barrier_epsilon);

        // Convert from the boolean vector to a vector of free dof indices
        Eigen::VectorXi free_dof(
            barrier_problem.num_vars - barrier_problem.fixed_dof.count());
        for (int i = 0, j = 0; i < barrier_problem.fixed_dof.size(); i++) {
            if (!barrier_problem.fixed_dof(i)) {
                free_dof(j++) = i;
            }
        }

        // Calculate the maximum number of iteration allowable
        int max_inner_iterations = int(max_iterations
            / ceil(-log2(min_barrier_epsilon) + log2(barrier_problem.epsilon)));

        barrier_problem.eval_intermediate_callback(barrier_problem.x0);

        OptimizationResults results;
        do {
            // Log the epsilon and the newton method will log the number of
            // iterations.
            spdlog::trace("solver=barrier_newton ϵ={:g}",
                barrier_constraint->barrier_epsilon);

            // Optimize for a fixed epsilon
            results
                = newtons_method(barrier_problem, free_dof, absolute_tolerance,
                    line_search_tolerance, max_inner_iterations);
            // Save the original problems objective
            results.minf = general_problem.eval_f(results.x);

            // Steepen the barrier
            barrier_problem.epsilon /= 2;
            barrier_constraint->barrier_epsilon = barrier_problem.epsilon;

            // Start next iteration from the ending optimal position
            barrier_problem.x0 = results.x;
        } while (barrier_problem.epsilon > min_barrier_epsilon);

        // TODO: This should check if the barrier constraints are satisfied.
        results.success = results.minf >= 0
            && barrier_problem.eval_f(results.x)
                < std::numeric_limits<double>::infinity();

        return results;
    }

} // namespace opt
} // namespace ccd
