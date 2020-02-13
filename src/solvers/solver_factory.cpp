#include "solver_factory.hpp"

#include <solvers/barrier_solver.hpp>
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
        barrier_inner_solvers_.emplace(
            "newton_solver", std::make_shared<NewtonSolver>("newton_solver"));
        barrier_inner_solvers_.emplace(
            "gradient_descent_solver",
            std::make_shared<GradientDescentSolver>("gradient_descent_solver"));
    }

    std::shared_ptr<IBarrierOptimizationSolver>
    SolverFactory::get_barrier_inner_solver(const std::string& problem) const
    {
        auto it = barrier_inner_solvers_.find(problem);
        assert(it != barrier_inner_solvers_.end());
        return it->second;
    }

} // namespace opt
} // namespace ccd
