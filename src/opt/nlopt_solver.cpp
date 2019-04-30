// Solve a optimization problem with NLopt.
#include <opt/nlopt_solver.hpp>
// #include <opt/solver.hpp>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>
#include <ccd/not_implemented_error.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    // Optimize the displacments using NLopt.
    OptimizationResults solve_problem_with_nlopt(
        OptimizationProblem& problem, const SolverSettings& settings)
    {
        nlopt::algorithm algorithm;
        switch (settings.method) {
        case MMA:
            algorithm = nlopt::LD_MMA;
            break;
        case SLSQP:
            algorithm = nlopt::LD_SLSQP;
            break;
        default:
            throw NotImplementedError(
                "Optimization method not implemented in NLopt.");
        }
        // Create a NLopt optimization object
        nlopt::opt opt(algorithm, unsigned(problem.num_vars));

        // Set the objective function
        opt.set_min_objective(nlopt_objective, &problem);

        // Set bound constraints
        opt.set_lower_bounds(std::vector<double>(problem.x_lower.data(),
            problem.x_lower.data() + problem.x_lower.size()));
        opt.set_upper_bounds(std::vector<double>(problem.x_upper.data(),
            problem.x_upper.data() + problem.x_upper.size()));

        // Set inequality constraints if desired
        std::vector<double> tol(2 * problem.num_constraints,
            settings.absolute_tolerance); // Tolerances
        // Set the m-constraints function
        opt.add_inequality_mconstraint(
            nlopt_inequality_constraints, &problem, tol);

        // Stopping criteria
        opt.set_xtol_rel(settings.relative_tolerance);
        // This may not play well with MMA for some reason
        opt.set_xtol_abs(settings.absolute_tolerance);
        opt.set_maxeval(settings.max_iter);
        opt.set_maxtime(settings.max_time);

        // Initial guess is the value of Uopt
        std::vector<double> x0(
            problem.x0.data(), problem.x0.data() + problem.x0.size());
        double minf; // The minimum objective value, upon return

        // Optimize the displacments
        nlopt::result r = opt.optimize(x0, minf);
        if (settings.verbosity) {
            print_nlopt_termination_reason(opt, r);
        }

        OptimizationResults results;
        results.x = Eigen::Map<Eigen::VectorXd>(x0.data(), problem.num_vars);
        results.minf = minf;
        Eigen::ArrayXd gx = problem.g(results.x).array();
        results.success = r > 0 && minf >= 0
            && problem.are_constraints_satisfied(
                results.x, settings.absolute_tolerance);
        return results;
    }

    // Computes NLopt's objective for an OptimizationProblem.
    double nlopt_objective(
        const std::vector<double>& x, std::vector<double>& grad, void* data)
    {
        assert(data); // Need the problem to compute anything.
        OptimizationProblem* problem = static_cast<OptimizationProblem*>(data);

        // Map the input to an Eigen Vector
        const Eigen::VectorXd X
            = Eigen::Map<const Eigen::VectorXd>(x.data(), problem->num_vars);

        // Only compute the gradient as needed.
        if (!grad.empty()) {
            // This should always be true according to NLopt.
            assert(size_t(grad.size()) == size_t(problem->num_vars));

            // Store the gradient of f(X) in grad
            Eigen::VectorXd::Map(grad.data(), grad.size()) = problem->grad_f(X);

            // This should really be done in the outer loop of NLopt (as Daniele
            // pointed out), but that would require modifying NLopt's source
            // code.
            // problem->intermediate_cb();
        }

        // Compute the objective at x
        return problem->f(X);
    }

    // Computes NLopt's constraints for an OptimizationProblem
    void nlopt_inequality_constraints(unsigned m, double* results, unsigned n,
        const double* x, double* grad, void* data)
    {
        assert(data); // Need data to compute volumes
        OptimizationProblem* problem = static_cast<OptimizationProblem*>(data);

        // Map the input to an Eigen Vector
        const Eigen::MatrixXd X = Eigen::Map<const Eigen::VectorXd>(x, n);

        // Compute the constraints at x
        Eigen::VectorXd gx = problem->g(X);

        // We want g_lower ≤ g(x) ≤ g_upper, but NLopt expects all
        // constraints(x) ≤ 0, so stack the constraints and negate as needed.
        assert(2 * problem->num_constraints == m);
        Eigen::VectorXd::Map(results, problem->num_constraints)
            = -gx + problem->g_lower;
        Eigen::VectorXd::Map(
            results + problem->num_constraints, problem->num_constraints)
            = gx - problem->g_upper;

        // Only compute the gradient as needed.
        if (grad) {
            // Gradient of the volumes (n x m)
            Eigen::MatrixXd dgx = problem->jac_g(X).transpose();
            // Check that the sizes match up
            assert(size_t(m * n) == size_t(2 * dgx.size()));
            // Flatten the Jacobian
            dgx.resize(dgx.size(), 1);
            // Store the flattened Jacobian in grad (once for each bound)
            Eigen::VectorXd::Map(grad, dgx.size()) = -dgx;
            Eigen::VectorXd::Map(grad + dgx.size(), dgx.size()) = dgx;

            // This should really be done in the outer loop of NLopt (as Daniele
            // pointed out), but that would require modifying NLopt's source
            // code.
            // problem->intermediate_cb();
        }
    } // namespace opt

    // Output the reason for NLopt's termination.
    void print_nlopt_termination_reason(
        const nlopt::opt& opt, const nlopt::result result)
    {
        switch (result) {
        case nlopt::XTOL_REACHED:
            std::cout << "Optimization terminated because the relative change "
                         "was less than "
                      << opt.get_xtol_rel() << "." << std::endl;
            break;
        case nlopt::MAXEVAL_REACHED:
            std::cout << "Optimization terminated because the number of "
                         "evaluations exceeded "
                      << opt.get_maxeval() << "." << std::endl;
            break;
        case nlopt::MAXTIME_REACHED:
            std::cout << "Optimization terminated because the runtime exceeded "
                      << opt.get_maxtime() << " seconds." << std::endl;
            break;
        default:
            std::cout << "Optimization terminated because of code " << result
                      << " (see "
                         "https://nlopt.readthedocs.io/en/latest/"
                         "NLopt_Reference/#return-values)."
                      << std::endl;
            break;
        }
    }

} // namespace opt
} // namespace ccd
