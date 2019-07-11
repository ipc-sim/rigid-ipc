#pragma once

#include <memory> // shared_ptr

#include <physics/simulation_problem.hpp>

namespace ccd {
namespace physics {
    class ProblemFactory {
    public:
        static const ProblemFactory& factory();

        std::shared_ptr<SimulationProblem> get_problem(const std::string& problem) const;
        inline const std::vector<std::string>& get_problem_names() const
        {
            return problem_names_;
        }

    private:
        ProblemFactory();
        std::map<std::string, std::shared_ptr<SimulationProblem>> problems_;
        std::vector<std::string> problem_names_;
    };
} // namespace opt

} // namespace ccd
