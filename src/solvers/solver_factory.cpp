#include "solver_factory.hpp"

#include <solvers/barrier_solver.hpp>
#include <solvers/bfgs_solver.hpp>
#include <solvers/gradient_descent_solver.hpp>
#include <solvers/ncp_solver.hpp>
#include <solvers/newton_solver.hpp>

namespace ccd {
namespace opt {

    const SolverFactory& SolverFactory::factory()
    {
        static SolverFactory instance;

        return instance;
    }

    SolverFactory::SolverFactory()
    {
        solvers_.emplace("barrier_solver",
            std::make_shared<BarrierSolver>("barrier_solver"));
        solvers_.emplace(
            "ncp_solver", std::make_shared<NCPSolver>("ncp_solver"));
        // inner solvers
        solvers_.emplace(
            "newton_solver", std::make_shared<NewtonSolver>("newton_solver"));
        solvers_.emplace(
            "bfgs_solver", std::make_shared<BFGSSolver>("bfgs_solver"));
        solvers_.emplace("gradient_descent_solver",
            std::make_shared<GradientDescentSolver>("gradient_descent_solver"));

        for (auto it = solvers_.begin(); it != solvers_.end(); ++it) {
            solver_names_.push_back(it->first);
        }
    }

    std::shared_ptr<OptimizationSolver> SolverFactory::get_solver(
        const std::string& problem) const
    {
        auto it = solvers_.find(problem);
        assert(it != solvers_.end());
        return it->second;
    }

} // namespace opt
} // namespace ccd
