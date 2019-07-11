#include "solver_factory.hpp"

#include <solvers/barrier_solver.hpp>
#include <solvers/ncp_solver.hpp>

namespace ccd {
namespace opt {

    const SolverFactory& SolverFactory::factory()
    {
        static SolverFactory instance;

        return instance;
    }

    SolverFactory::SolverFactory()
    {
        solvers_.emplace("barrier_solver", std::make_shared<BarrierSolver>());
        solvers_.emplace("ncp_solver", std::make_shared<NCPSolver>());

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
