#pragma once
#include <memory> // shared_ptr

#include <solvers/barrier_solver.hpp>

namespace ccd {
namespace opt {
    class SolverFactory {
    public:
        static const SolverFactory& factory();

        std::shared_ptr<BarrierInnerSolver>
        get_barrier_inner_solver(const std::string& problem) const;

    private:
        SolverFactory();

        std::map<std::string, std::shared_ptr<BarrierInnerSolver>>
            barrier_inner_solvers;
    };
} // namespace opt
} // namespace ccd
