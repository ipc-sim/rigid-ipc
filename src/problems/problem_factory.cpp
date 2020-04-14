#include "problem_factory.hpp"

#include <solvers/barrier_solver.hpp>

#include <problems/distance_barrier_rb_problem.hpp>
// #include <problems/volume_rb_problem.hpp>

namespace ccd {

const ProblemFactory& ProblemFactory::factory()
{
    static ProblemFactory instance;

    return instance;
}

ProblemFactory::ProblemFactory()
{
    problems_.emplace(
        opt::DistanceBarrierRBProblem::problem_name(),
        std::make_shared<opt::DistanceBarrierRBProblem>());
    // problems_.emplace(
    //     "volume_rb_problem", std::make_shared<opt::VolumeRBProblem>());
}

std::shared_ptr<physics::SimulationProblem>
ProblemFactory::get_problem(const std::string& name) const
{
    auto it = problems_.find(name);
    assert(it != problems_.end());
    return it->second;
}

} // namespace ccd
