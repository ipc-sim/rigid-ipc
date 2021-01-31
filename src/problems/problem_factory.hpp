#pragma once

#include <memory> // shared_ptr

#include <physics/simulation_problem.hpp>
#include <solvers/optimization_solver.hpp>

namespace ipc::rigid {

class ProblemFactory {
public:
    static const ProblemFactory& factory();

    std::shared_ptr<SimulationProblem>
    get_problem(const std::string& name) const;

private:
    ProblemFactory();
    std::map<std::string, std::shared_ptr<SimulationProblem>>
        problems_;
};

} // namespace ipc::rigid
