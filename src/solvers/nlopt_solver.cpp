// Solve a optimization problem with NLopt.
#include "nlopt_solver.hpp"

#include <utils/not_implemented_error.hpp>

#include <iostream>

namespace ccd {
namespace opt {

    NLOptSolver::NLOptSolver()
        : algorithm(nlopt::LD_SLSQP)
        , absolute_tolerance(1e-8)
        , relative_tolerance(1e-8)
        , max_iterations(1000)
        , max_time(2e19)
        , verbose(false)
    {
    }

    NLOptSolver::~NLOptSolver() {}

    OptimizationResults NLOptSolver::solve(IConstraintedProblem& problem)
    {
        nlopt::opt opt(algorithm, unsigned(problem.num_vars()));

        opt.set_min_objective(nlopt_objective, &problem);

        // Set inequality constraints if desired
        std::vector<double> tol(
            uint(2 * problem.num_constraints()), absolute_tolerance);
        opt.add_inequality_mconstraint(
            nlopt_inequality_constraints, &problem, tol);

        // Stopping criteria
        opt.set_xtol_rel(relative_tolerance);
        // This may not play well with MMA for some reason
        opt.set_xtol_abs(absolute_tolerance);
        opt.set_maxeval(max_iterations);
        opt.set_maxtime(max_time);

        // Initial guess is the value of Uopt
        auto const& starting_point = problem.starting_point();
        std::vector<double> x0(starting_point.data(),
            starting_point.data() + starting_point.size());
        double minf; // The minimum objective value, upon return

        // Optimize the displacements
        nlopt::result r = opt.optimize(x0, minf);
        if (verbose) {
            print_nlopt_termination_reason(opt, r);
        }

        OptimizationResults results;
        results.x = Eigen::Map<Eigen::VectorXd>(x0.data(), problem.num_vars());
        results.minf = minf;
        Eigen::ArrayXd gx = problem.eval_g(results.x).array();
        results.success = r > 0 && minf >= 0;
        return results;
    }

    double nlopt_objective(
        const std::vector<double>& x, std::vector<double>& grad, void* data)
    {
        assert(data != nullptr);
        IConstraintedProblem* problem = static_cast<IConstraintedProblem*>(data);

        const Eigen::VectorXd X
            = Eigen::Map<const Eigen::VectorXd>(x.data(), problem->num_vars());

        if (!grad.empty()) {
            // This should always be true according to NLopt.
            assert(size_t(grad.size()) == size_t(problem->num_vars()));

            // Store the gradient of f(X) in grad
            Eigen::VectorXd::Map(grad.data(), int(grad.size()))
                = problem->eval_grad_f(X);
        }

        return problem->eval_f(X);
    }

    void nlopt_inequality_constraints(unsigned m,
        double* results,
        unsigned n,
        const double* x,
        double* grad,
        void* data)
    {
        assert(data); // Need data to compute volumes
        IConstraintedProblem* problem = static_cast<IConstraintedProblem*>(data);

        const Eigen::MatrixXd X = Eigen::Map<const Eigen::VectorXd>(x, n);
        Eigen::VectorXd gx = problem->eval_g(X);

        // We want g(x) >= 0, but NLopt expects all
        // constraints(x) ≤ 0, so stack the constraints and negate as needed.
        assert(uint(2 * problem->num_constraints()) == m);
        Eigen::VectorXd::Map(results, problem->num_constraints()) = -gx;

        if (grad) {
            Eigen::MatrixXd dgx = problem->eval_jac_g(X).transpose();
            assert(size_t(m * n) == size_t(2 * dgx.size()));
            // Flatten the Jacobian
            dgx.resize(dgx.size(), 1);
            // Store the flattened Jacobian in grad (once for each bound)
            Eigen::VectorXd::Map(grad, dgx.size()) = -dgx;
            Eigen::VectorXd::Map(grad + dgx.size(), dgx.size()) = dgx;
        }
    }

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
