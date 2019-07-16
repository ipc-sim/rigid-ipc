#include "problem_factory.hpp"

#include <physics/particles_problem.hpp>
#include <physics/rigid_body_problem.hpp>

namespace ccd {
namespace physics {

    const ProblemFactory& ProblemFactory::factory()
    {
        static ProblemFactory instance;

        return instance;
    }

    ProblemFactory::ProblemFactory()
    {
        problems_.emplace("rigid_body_problem",
            std::make_shared<RigidBodyProblem>("rigid_body_problem"));
        problems_.emplace("particles_problem",
            std::make_shared<ParticlesDisplProblem>("particles_problem"));

        for (auto it = problems_.begin(); it != problems_.end(); ++it) {
            problem_names_.push_back(it->first);
        }
    }

    std::shared_ptr<SimulationProblem> ProblemFactory::get_problem(
        const std::string& problem) const
    {
        auto it = problems_.find(problem);
        assert(it != problems_.end());
        return it->second;
    }

} // namespace physics
} // namespace ccd
