#include "problem_factory.hpp"

#include <solvers/barrier_solver.hpp>

#include <physics/particles_problem.hpp>
#include <problems/distance_barrier_rb_problem.hpp>

namespace ccd {

const ProblemFactory& ProblemFactory::factory()
{
    static ProblemFactory instance;

    return instance;
}

ProblemFactory::ProblemFactory()
{
    problems_.emplace("distance_barrier_rb_problem",
        std::make_shared<opt::DistanceBarrierRBProblem>(
            "distance_barrier_rb_problem"));
}

std::shared_ptr<physics::ISimulationProblem> ProblemFactory::get_problem(
    const std::string& name) const
{
    auto it = problems_.find(name);
    assert(it != problems_.end());
    return it->second;
}

} // namespace ccd
