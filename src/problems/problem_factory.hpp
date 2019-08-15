#pragma once

#include <memory> // shared_ptr

#include <physics/simulation_problem.hpp>
#include <solvers/optimization_solver.hpp>

namespace ccd {

class ProblemFactory {
public:
    static const ProblemFactory& factory();

    std::shared_ptr<physics::ISimulationProblem> get_problem(
        const std::string& name) const;

private:
    ProblemFactory();
    std::map<std::string, std::shared_ptr<physics::ISimulationProblem>>
        problems_;
};

} // namespace ccd
