#include "solver_factory.hpp"

#include <solvers/barrier_solver.hpp>
#include <solvers/gradient_descent_solver.hpp>
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
        barrier_inner_solvers.emplace(
            NewtonSolver::solver_name(), std::make_shared<NewtonSolver>());
        barrier_inner_solvers.emplace(
            GradientDescentSolver::solver_name(),
            std::make_shared<GradientDescentSolver>());
    }

    std::shared_ptr<BarrierInnerSolver>
    SolverFactory::get_barrier_inner_solver(const std::string& problem) const
    {
        auto it = barrier_inner_solvers.find(problem);
        assert(it != barrier_inner_solvers.end());
        return it->second;
    }

} // namespace opt
} // namespace ccd
