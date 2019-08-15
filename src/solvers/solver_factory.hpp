#pragma once
#include <memory> // shared_ptr

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {
    class SolverFactory {
    public:
        static const SolverFactory& factory();

        std::shared_ptr<IBarrierOptimizationSolver> get_barrier_inner_solver(
            const std::string& problem) const;


    private:
        SolverFactory();

        std::map<std::string, std::shared_ptr<IBarrierOptimizationSolver>>
            barrier_inner_solvers_;

    };
} // namespace opt
} // namespace ccd
