// Methods for optimizing the displacments with a non-linear interference volume
// constraint.
#include <opt/displacement_opt.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/ncp_solver.hpp>

#include <logger.hpp>

namespace ccd {
namespace opt {

    // Optimize the displacment opt problem with the given method and starting
    // value.
    OptimizationResults displacement_optimization(OptimizationProblem& problem,
        const Eigen::MatrixX2d& U0, SolverSettings& settings)
    {

        // initial value
        Eigen::MatrixXd x0 = U0;
        x0.resize(U0.size(), 1); // Flatten displacements
        problem.x0 = x0;

        OptimizationResults result;
        if (settings.method == NCP) {
            result = solve_ncp_displacement_optimization(problem, settings);
        } else {
            result = solve_problem(problem, settings);
        }
        result.x.resize(U0.rows(), 2); // Unflatten displacments
        return result;
    } // namespace opt

    OptimizationResults solve_ncp_displacement_optimization(
        OptimizationProblem& problem, SolverSettings& settings)
    {
        // Solves the KKT conditions of the Optimization Problem
        //  (U - Uk) = \nabla g(U)
        //  s.t V(U) >= 0
        int num_vars = int(problem.x0.rows());
        Eigen::SparseMatrix<double> A(num_vars, num_vars);
        A.setIdentity();

        Eigen::VectorXd b, x_opt, lambda_opt;
        // obtain Uk from grad_f using x=0
        b = -problem.grad_f(Eigen::VectorXd::Zero(num_vars));

        int num_it = 0;
        callback_intermediate_ncp callback =
            [&problem, &settings, &num_it](const Eigen::VectorXd& x,
                const Eigen::VectorXd& alpha, const double gamma) {
                settings.intermediate_cb(x, problem.f(x), alpha, gamma, num_it);
                num_it += 1;
            };
        OptimizationResults result;
        result.success
            = solve_ncp(A, b, problem.g, problem.jac_g, settings.max_iter,
                callback, NcpUpdate::G_GRADIENT, settings.lcp_solver, x_opt,
                lambda_opt, /*check_convergence=*/false);
        result.x = x_opt;
        result.minf = problem.f(x_opt);

        return result;
    }

} // namespace opt
} // namespace ccd
