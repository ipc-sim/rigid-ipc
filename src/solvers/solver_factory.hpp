#pragma once
#include <memory> // shared_ptr

#include <solvers/optimization_solver.hpp>

namespace ipc::rigid {

class SolverFactory {
public:
    static const SolverFactory& factory();

    std::shared_ptr<OptimizationSolver>
    get_barrier_solver(const std::string& solver_name) const;

private:
    SolverFactory();

    std::map<std::string, std::shared_ptr<OptimizationSolver>> barrier_solvers;
};

} // namespace ipc::rigid
