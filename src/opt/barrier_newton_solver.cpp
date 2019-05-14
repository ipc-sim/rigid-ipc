// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include <opt/barrier_newton_solver.hpp>

#include <iostream>

#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

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
        auto eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier(xi, eps); });
    }

    Eigen::VectorXd BarrierProblem::barrier_gradient(const Eigen::VectorXd x)
    {
        auto eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier_gradient(xi, eps); });
    }

    Eigen::VectorXd BarrierProblem::barrier_hessian(const Eigen::VectorXd x)
    {
        auto eps = epsilon;
        return x.unaryExpr(
            [eps](double xi) { return spline_barrier_hessian(xi, eps); });
    }
    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {

        Eigen::VectorXd gx = general_problem->eval_g(x);
        double val = general_problem->eval_f(x)
            + barrier(x - general_problem->x_lower).sum()
            + barrier(-x + general_problem->x_upper).sum();

        // Check to make sure the gx is not dynamic
        if (gx.size() != 0
            && gx.size() == gx.rows() == general_problem->g_lower.rows()
            && gx.cols() == general_problem->g_lower.cols()
            && gx.rows() == general_problem->g_upper.rows()
            && gx.cols() == general_problem->g_upper.cols()) {
            val += barrier(gx - general_problem->g_lower).sum()
                + barrier(-gx + general_problem->g_upper).sum();
        } else {
            val += barrier(gx).sum();
        }

        return val;
    }

    Eigen::VectorXd BarrierProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd gx = general_problem->eval_g(x);
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);
        Eigen::VectorXd grad = general_problem->eval_grad_f(x);

        // Add constraint functions barrier(g(x))
        // Check to make sure the gx is not dynamic
        Eigen::VectorXd coeffs;
        if (gx.size() != 0 && dgx.size() != 0
            && gx.rows() == general_problem->g_lower.rows()
            && gx.cols() == general_problem->g_lower.cols()
            && gx.rows() == general_problem->g_upper.rows()
            && gx.cols() == general_problem->g_upper.cols()) {
            coeffs = barrier_gradient(gx - general_problem->g_lower)
                - barrier_gradient(-gx + general_problem->g_upper);
        } else {
            coeffs = barrier_gradient(gx);
        }
        grad += (coeffs.asDiagonal() * dgx).colwise().sum().transpose();

        // Add value constraints barrier(x)
        // âˆ‡ âˆ‘ Ï•(x_i) = âˆ‘ (âˆ‡ Ï•(x_i)) = âˆ‘ [0 ... Ï•'(x_i) ... 0]^T
        //            = [Ï•'(x_1) Ï•'(x_2) ... Ï•'(x_n)]^T
        grad += barrier_gradient(x - general_problem->x_lower)
            - barrier_gradient(-x + general_problem->x_upper);

        return grad;
    }

    Eigen::MatrixXd BarrierProblem::eval_hessian_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd gx = general_problem->eval_g(x);
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);
        std::vector<Eigen::MatrixXd> ddgx = general_problem->eval_hessian_g(x);

        Eigen::MatrixXd hessian = general_problem->eval_hessian_f(x);

        // Check to make sure the gx is not dynamic
        Eigen::VectorXd grad_coeffs, hessian_coeffs;
        if (gx.size() != 0 && dgx.size() != 0 && ddgx.size() != 0
            && gx.rows() == general_problem->g_lower.rows()
            && gx.cols() == general_problem->g_lower.cols()
            && gx.rows() == general_problem->g_upper.rows()
            && gx.cols() == general_problem->g_upper.cols()) {
            grad_coeffs = barrier_gradient(gx - general_problem->g_lower)
                - barrier_gradient(-gx + general_problem->g_upper);
            hessian_coeffs = barrier_hessian(gx - general_problem->g_lower)
                + barrier_hessian(-gx + general_problem->g_upper);
        } else {
            grad_coeffs = barrier_gradient(gx);
            hessian_coeffs = barrier_hessian(gx);
        }

        for (long i = 0; i < gx.rows(); i++) {
            hessian += hessian_coeffs(i) * dgx.row(i).transpose() * dgx.row(i)
                + grad_coeffs(i) * ddgx[unsigned(i)];
        }

        // âˆ‡ [Ï•'(x_1) Ï•'(x_2) ... Ï•'(x_n)]^T
        // = diag([Ï•''(x_1) Ï•''(x_2) ... Ï•''(x_n)]^T)
        hessian.diagonal() += barrier_hessian(x - general_problem->x_lower)
            + barrier_hessian(-x + general_problem->x_upper);
        return hessian;
    }

    BarrierNewtonSolver::BarrierNewtonSolver()
        : barrier_epsilon(0.0)
        , min_barrier_epsilon(1e-8)
        , absolute_tolerance(1e-8)
        , line_search_tolerance(1e-8)
        , max_iterations(3000)
    {
    }

    BarrierNewtonSolver::~BarrierNewtonSolver() {}

    OptimizationResults BarrierNewtonSolver::solve(
        OptimizationProblem& general_problem)
    {
        BarrierProblem barrier_problem(general_problem, barrier_epsilon);

        OptimizationResults results;
        do {
            results = newtons_method(barrier_problem, absolute_tolerance,
                line_search_tolerance, max_iterations);
            // Save the original problems objective
            results.minf = general_problem.eval_f(results.x);
            // Steepen the barrier
            barrier_problem.epsilon /= 2;
            // Start next iteration from the ending optimal position
            barrier_problem.x0 = results.x;
        } while (barrier_problem.epsilon > min_barrier_epsilon);

        results.success = results.minf >= 0
            && general_problem.are_constraints_satisfied(results.x, 0.0);

        return results;
    }


} // namespace opt
} // namespace ccd
