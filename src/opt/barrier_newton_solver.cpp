// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include <opt/barrier_newton_solver.hpp>

#include <iostream>

#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

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

        f_uk += gx.sum() + this->barrier(x - general_problem->x_lower).sum()
            + this->barrier(-x + general_problem->x_upper).sum();
        f_uk_gradient += dgx.colwise().sum().transpose()
            + barrier_gradient(x - general_problem->x_lower)
            - barrier_gradient(-x + general_problem->x_upper);

        PROFILE(
            for (const auto& ddgx_i
                 : ddgx) { f_uk_hessian += ddgx_i; },
            ProfiledPoint::SUMMING_HESSIAN)
    }

    double BarrierProblem::eval_f(const Eigen::VectorXd& x)
    {
        double val = general_problem->eval_f(x)
            + this->barrier(x - general_problem->x_lower).sum()
            + this->barrier(-x + general_problem->x_upper).sum();

        Eigen::VectorXd gx = general_problem->eval_g(x);
        // Check to make sure the gx is not dynamic
        // if (gx.size() != 0
        //     && gx.size() == gx.rows() == general_problem->g_lower.rows()
        //     && gx.cols() == general_problem->g_lower.cols()
        //     && gx.rows() == general_problem->g_upper.rows()
        //     && gx.cols() == general_problem->g_upper.cols()) {
        //     val += barrier(gx - general_problem->g_lower).sum()
        //         + barrier(-gx + general_problem->g_upper).sum();
        // } else {
        //     val += barrier(gx).sum();
        // }
        val += gx.sum();

        return val;
    }

    Eigen::VectorXd BarrierProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd grad = general_problem->eval_grad_f(x);

        // Add constraint functions barrier(g(x))
        // Eigen::VectorXd gx = general_problem->eval_g(x);
        Eigen::MatrixXd dgx = general_problem->eval_jac_g(x);
        // Check to make sure the gx is not dynamic
        // Eigen::VectorXd coeffs;
        // if (gx.size() != 0 && dgx.size() != 0
        //     && gx.rows() == general_problem->g_lower.rows()
        //     && gx.cols() == general_problem->g_lower.cols()
        //     && gx.rows() == general_problem->g_upper.rows()
        //     && gx.cols() == general_problem->g_upper.cols()) {
        //     coeffs = barrier_gradient(gx - general_problem->g_lower)
        //         - barrier_gradient(-gx + general_problem->g_upper);
        // } else {
        //     coeffs = barrier_gradient(gx);
        // }
        // grad += (coeffs.asDiagonal() * dgx).colwise().sum().transpose();
        grad += dgx.colwise().sum().transpose();

        // Add value constraints barrier(x)
        // ∇ ∑ ϕ(x_i) = ∑ (∇ ϕ(x_i)) = ∑ [0 ... ϕ'(x_i) ... 0]^T
        //            = [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
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

        // Check to make sure the gx is not dynamic
        // Eigen::VectorXd grad_coeffs, hessian_coeffs;
        // if (gx.size() != 0 && dgx.size() != 0 && ddgx.size() != 0
        //     && gx.rows() == general_problem->g_lower.rows()
        //     && gx.cols() == general_problem->g_lower.cols()
        //     && gx.rows() == general_problem->g_upper.rows()
        //     && gx.cols() == general_problem->g_upper.cols()) {
        //     grad_coeffs = barrier_gradient(gx - general_problem->g_lower)
        //         - barrier_gradient(-gx + general_problem->g_upper);
        //     hessian_coeffs = barrier_hessian(gx - general_problem->g_lower)
        //         + barrier_hessian(-gx + general_problem->g_upper);
        // } else {
        //     grad_coeffs = barrier_gradient(gx);
        //     hessian_coeffs = barrier_hessian(gx);
        // }

        PROFILE(
            if (ddgx.size() > 0) {
                Eigen::SparseMatrix<double> sum_ddgx = ddgx[0];
                for (long i = 1; i < ddgx.size(); i++) {
                    // hessian += hessian_coeffs(i) * dgx.row(i).transpose() *
                    // dgx.row(i)
                    //     + grad_coeffs(i) * ddgx[unsigned(i)];
                    sum_ddgx += ddgx[unsigned(i)];
                }
                hessian += sum_ddgx;
            },
            ProfiledPoint::SUMMING_HESSIAN)

        // ∇ [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
        // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
        hessian.diagonal() += barrier_hessian(x - general_problem->x_lower)
            + barrier_hessian(-x + general_problem->x_upper);
        return hessian;
    }

    Eigen::SparseMatrix<double> BarrierProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        Eigen::SparseMatrix<double> hessian
            = general_problem->eval_hessian_f_sparse(x);

        std::vector<Eigen::SparseMatrix<double>> ddgx
            = general_problem->eval_hessian_g(x);

        PROFILE(
            for (const auto& ddgx_i
                 : ddgx) { hessian += ddgx_i; },
            ProfiledPoint::SUMMING_HESSIAN)

        // ∇ [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
        // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
        // hessian.diagonal() += barrier_hessian(x - general_problem->x_lower)
        //     + barrier_hessian(-x + general_problem->x_upper);
        return hessian;
    }

    bool BarrierProblem::eval_intermediate_callback(const Eigen::VectorXd& x)
    {
        return general_problem->eval_intermediate_callback(x);
    }

    BarrierNewtonSolver::BarrierNewtonSolver()
        : barrier_constraint(nullptr)
        , min_barrier_epsilon(1e-8)
        , absolute_tolerance(1e-8)
        , line_search_tolerance(1e-8)
        , max_iterations(3000)
        , verbose(false)
    {
    }

    BarrierNewtonSolver::~BarrierNewtonSolver() {}

    OptimizationResults BarrierNewtonSolver::solve(
        OptimizationProblem& general_problem)
    {
        assert(barrier_constraint != nullptr);

        BarrierProblem barrier_problem(
            general_problem, barrier_constraint->barrier_epsilon);

        if (verbose) {
            std::cout << std::endl;
        }

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

        OptimizationResults results;
        while (barrier_problem.epsilon > min_barrier_epsilon) {
            // Log the epsilon and the newton method will log the number of
            // iterations.
            if (verbose) {
                std::cout << "ϵ = " << barrier_constraint->barrier_epsilon
                          << ": ";
            }

            // Optimize for a fixed epsilon
            results
                = newtons_method(barrier_problem, free_dof, absolute_tolerance,
                    line_search_tolerance, max_inner_iterations, verbose);
            // Save the original problems objective
            results.minf = general_problem.eval_f(results.x);

            // Steepen the barrier
            barrier_problem.epsilon /= 2;
            barrier_constraint->barrier_epsilon = barrier_problem.epsilon;

            // Start next iteration from the ending optimal position
            barrier_problem.x0 = results.x;
        }

        results.success = results.minf >= 0
            && general_problem.are_constraints_satisfied(results.x, 0.0);

        return results;
    }

} // namespace opt
} // namespace ccd
