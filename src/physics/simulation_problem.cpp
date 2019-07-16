#include "simulation_problem.hpp"

namespace ccd {
namespace physics {
    SimulationProblem::SimulationProblem(const std::string& name)
        : opt::OptimizationProblem(name)
    {
    }

} // namespace physics
} // namespace ccd
