// Functions for optimizing functions.
#include "bfgs_solver.hpp"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/slice.h>
#include <igl/slice_into.h>

#include <solvers/lbfgspp/LBFGS.h>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    BFGSSolver::BFGSSolver()
        : BFGSSolver("bfgs_solver")
    {
    }

    BFGSSolver::BFGSSolver(const std::string& name)
        : absolute_tolerance(1e-5)
        , min_step_length(1e-12)
        , name_(name)
    {
    }

    BFGSSolver::~BFGSSolver() {}

    void BFGSSolver::settings(const nlohmann::json& json)
    {
        max_iterations = json["max_iterations"].get<int>();
        absolute_tolerance = json["absolute_tolerance"].get<double>();
        min_step_length = json["min_step_length"].get<double>();
    }

    nlohmann::json BFGSSolver::settings() const
    {
        nlohmann::json json;
        json["max_iterations"] = max_iterations;
        json["absolute_tolerance"] = absolute_tolerance;
        json["min_step_length"] = min_step_length;
        return json;
    }

    OptimizationResults BFGSSolver::solve(IBarrierProblem& problem)
    {
        // Functional in the form expected by LBFGS++
        auto f = [&problem](const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
            double fx;
            problem.eval_f_and_fdiff(x, fx, grad);
            return fx;
        };

        // Callback as a function pointer
        auto callback = [&problem](const Eigen::VectorXd& x) -> bool {
            // TODO: refactor to allow intermediate steps instead
            return true;
        };

        // Set up parameters
        LBFGSpp::LBFGSParam<double> param;
        param.epsilon = absolute_tolerance;
        param.min_step = min_step_length;
        param.max_iterations = max_iterations;
        param.max_linesearch = 60;
        param.linesearch = LBFGSpp::LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;

        // Create solver and function object
        LBFGSpp::LBFGSSolver<double> solver(param);

        std::string exit_reason = "found a local optimum";
        Eigen::VectorXd x = problem.starting_point();
        int niter = -1;
        try {
            double _;
            niter = solver.minimize(f, x, _, callback);
        } catch (std::runtime_error err) {
            exit_reason = err.what();
        }
        if (niter >= max_iterations) {
            exit_reason = "exceeded the maximum allowable iterations";
        }

        Eigen::VectorXd grad;
        f(x, grad);
        spdlog::trace("solver=BFGS total_iter={:d} exit_reason=\"{:s}\" "
                      "sqr_norm_grad={:g}",
            niter, exit_reason, grad.squaredNorm());

        OptimizationResults results;
        results.x = x;
        results.minf = problem.eval_f(x);
        results.success = !std::isinf(results.minf);
        if (!results.success) {
            spdlog::warn("solver=BFGS failure=\"results violate constraint "
                         "(f(x) < ∞)\"");
        }
        return results;
    }

    void BFGSSolver::init_free_dof(Eigen::VectorXb is_dof_fixed)
    {
        free_dof = Eigen::VectorXi(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
    }

} // namespace opt
} // namespace ccd
