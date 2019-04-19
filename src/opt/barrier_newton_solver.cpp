// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include <opt/barrier_newton_solver.hpp>

#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

namespace ccd {
namespace opt {

    // Performa Newton's Method to minimize a function f(x) + barrier(g(x)).
    OptimizationResults solve_problem_with_barrier_newton(
        const OptimizationProblem& general_problem, SolverSettings& settings)
    {
        OptimizationProblem barrier_problem;
        setup_barrier_problem(
            general_problem, settings.barrier_epsilon, barrier_problem);

        OptimizationResults results;
        do {
            results = newtons_method(barrier_problem, settings);
            // Save the original problems objective
            results.minf = general_problem.f(results.x);
            // Steepen the barrier
            settings.barrier_epsilon /= 2;
            // Start next iteration from the ending optimal position
            barrier_problem.x0 = results.x;
        } while (settings.barrier_epsilon > settings.min_barrier_epsilon);

        results.success = results.minf >= 0
            && general_problem.are_constraints_satisfied(results.x, 0.0);
        return results;
    }

    // Create a OptimizationProblem for the barrier Newton's method.
    void setup_barrier_problem(const OptimizationProblem& general_problem,
        double& epsilon, OptimizationProblem& barrier_problem)
    {
        barrier_problem.num_vars = general_problem.num_vars;
        barrier_problem.num_constraints = general_problem.num_constraints;
        barrier_problem.x0 = general_problem.x0;

        // Define the barrier function or a fixed epsilon
        auto barrier
            = [&epsilon](double x) { return spline_barrier(x, epsilon); };
        auto barrier_gradient = [&epsilon](double x) {
            return spline_barrier_gradient(x, epsilon);
        };
        auto barrier_hessian = [&epsilon](double x) {
            return spline_barrier_hessian(x, epsilon);
        };

        // Redefine the objective with the constraints as penalties
        barrier_problem.f = [&general_problem, barrier](
                                const Eigen::VectorXd& x) {
            Eigen::VectorXd gx = general_problem.g(x);
            double val = general_problem.f(x)
                + (x - general_problem.x_lower).unaryExpr(barrier).sum()
                + (-x + general_problem.x_upper).unaryExpr(barrier).sum();
            if (gx.size() != 0) {
                // Check to make sure the gx is not dynamic
                assert(gx.rows() == general_problem.g_lower.rows()
                    && gx.cols() == general_problem.g_lower.cols());
                assert(gx.rows() == general_problem.g_upper.rows()
                    && gx.cols() == general_problem.g_upper.cols());

                val += (gx - general_problem.g_lower).unaryExpr(barrier).sum()
                    + (-gx + general_problem.g_upper).unaryExpr(barrier).sum();
            }
            return val;
        };

        // Redefine the objective gradient with the constraints as penalties
        barrier_problem.grad_f = [&general_problem, barrier_gradient](
                                     const Eigen::VectorXd& x) {
            Eigen::VectorXd gx = general_problem.g(x);
            Eigen::MatrixXd dgx = general_problem.jac_g(x);

            Eigen::VectorXd grad = general_problem.grad_f(x);

            // Add constraint functions barrier(g(x))
            if (gx.size() != 0 && dgx.size() != 0) {
                // Check to make sure the gx is not dynamic
                assert(gx.rows() == general_problem.g_lower.rows()
                    && gx.cols() == general_problem.g_lower.cols());
                assert(gx.rows() == general_problem.g_upper.rows()
                    && gx.cols() == general_problem.g_upper.cols());

                Eigen::VectorXd coeffs
                    = (gx - general_problem.g_lower).unaryExpr(barrier_gradient)
                    - (-gx + general_problem.g_upper)
                          .unaryExpr(barrier_gradient);
                grad += (coeffs.asDiagonal() * dgx).colwise().sum().transpose();
            }

            // Add value constraints barrier(x)
            // ∇ ∑ ϕ(x_i) = ∑ (∇ ϕ(x_i)) = ∑ [0 ... ϕ'(x_i) ... 0]^T
            //            = [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
            grad += (x - general_problem.x_lower).unaryExpr(barrier_gradient)
                - (-x + general_problem.x_upper).unaryExpr(barrier_gradient);
            return grad;
        };

        // Redefine the objective hessian with the constraints as penalties
        barrier_problem.hessian_f = [&general_problem, barrier_gradient,
                                        barrier_hessian](
                                        const Eigen::VectorXd& x) {
            Eigen::VectorXd gx = general_problem.g(x);
            Eigen::MatrixXd dgx = general_problem.jac_g(x);
            std::vector<Eigen::MatrixXd> ddgx = general_problem.hessian_g(x);

            Eigen::MatrixXd hessian = general_problem.hessian_f(x);
            if (gx.size() != 0 && dgx.size() != 0 && ddgx.size() != 0) {
                // Check to make sure the gx is not dynamic
                assert(gx.rows() == general_problem.g_lower.rows()
                    && gx.cols() == general_problem.g_lower.cols());
                assert(gx.rows() == general_problem.g_upper.rows()
                    && gx.cols() == general_problem.g_upper.cols());

                Eigen::VectorXd hessian_coeffs
                    = (gx - general_problem.g_lower).unaryExpr(barrier_hessian)
                    + (-gx + general_problem.g_upper)
                          .unaryExpr(barrier_hessian);
                Eigen::VectorXd grad_coeffs
                    = (gx - general_problem.g_lower).unaryExpr(barrier_gradient)
                    - (-gx + general_problem.g_upper)
                          .unaryExpr(barrier_gradient);
                for (long i = 0; i < general_problem.num_constraints; i++) {
                    hessian += hessian_coeffs(i) * dgx.row(i).transpose()
                            * dgx.row(i)
                        + grad_coeffs(i) * ddgx[unsigned(i)];
                }
            }
            // ∇ [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
            // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
            hessian.diagonal()
                += (x - general_problem.x_lower).unaryExpr(barrier_hessian)
                + (-x + general_problem.x_upper).unaryExpr(barrier_hessian);
            return hessian;
        };
    }

} // namespace opt
} // namespace ccd
