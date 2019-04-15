// Solve the optimization general_problem using Newton's Method with barriers
// for the constraints.

#include <opt/barrier_newton_solver.hpp>

#include <opt/barrier.hpp>
#include <opt/newtons_method.hpp>

namespace ccd {
namespace opt {

    // Create a OptimizationProblem for displacment optimization
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
        barrier_problem.f
            = [&general_problem, barrier](const Eigen::VectorXd& x) {
                  auto gx = general_problem.g(x);
                  return general_problem.f(x)
                      + (gx - general_problem.g_lower).unaryExpr(barrier).sum()
                      + (-gx + general_problem.g_upper).unaryExpr(barrier).sum()
                      + (x - general_problem.x_lower).unaryExpr(barrier).sum()
                      + (-x + general_problem.x_upper).unaryExpr(barrier).sum();
              };

        // Redefine the objective gradient with the constraints as penalties
        barrier_problem.grad_f = [&general_problem, barrier_gradient](
                                     const Eigen::VectorXd& x) {
            auto gx = general_problem.g(x);
            auto dgx = general_problem.jac_g(x);

            auto grad = general_problem.grad_f(x);
            for (long i = 0; i < general_problem.num_constraints; i++) {
                grad += (barrier_gradient(gx(i) - general_problem.g_lower(i))
                            - barrier_gradient(
                                -gx(i) + general_problem.g_upper(i)))
                    * dgx.row(i).transpose();
            }
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
            auto gx = general_problem.g(x);
            auto dgx = general_problem.jac_g(x);
            auto ddgx = general_problem.hessian_g(x);

            auto hessian = general_problem.hessian_f(x);
            for (long i = 0; i < general_problem.num_constraints; i++) {
                hessian += (barrier_hessian(gx(i) - general_problem.g_lower(i))
                               + barrier_hessian(
                                   -gx(i) + general_problem.g_upper(i)))
                        * dgx.row(i).transpose() * dgx.row(i)
                    + (barrier_gradient(gx(i) - general_problem.g_lower(i))
                          - barrier_gradient(
                              -gx(i) + general_problem.g_upper(i)))
                        * ddgx[unsigned(i)];
            }
            // \nabla [ϕ'(x_1) ϕ'(x_2) ... ϕ'(x_n)]^T
            // = diag([ϕ''(x_1) ϕ''(x_2) ... ϕ''(x_n)]^T)
            hessian.diagonal()
                += (x - general_problem.x_lower).unaryExpr(barrier_hessian)
                + (-x + general_problem.x_upper).unaryExpr(barrier_hessian);
            return hessian;
        };
    }

    // Performa Newton's Method to minimize a function f(x).
    OptimizationResults solve_problem_with_barrier_newton(
        const OptimizationProblem& problem, const SolverSettings& settings)
    {
        // TODO: Replace this with maximal displacement over all vertices
        // double epsilon = general_problem.starting_epsilon;
        double epsilon = 1;

        // auto compute_objective[&](const Eigen::VectorXd& x,
        //     Eigen::VectorXd* gradient, Eigen::MatrixXd* hessian){}

        OptimizationProblem barrier_problem;
        setup_barrier_problem(problem, epsilon, barrier_problem);

        return newtons_method(barrier_problem, settings, /*mu = */ 1e-5);
    }

} // namespace opt
} // namespace ccd
