#pragma once
#include <memory> // shared_ptr

#include <solvers/optimization_solver.hpp>

namespace ccd {
namespace opt {
    class SolverFactory {
    public:
        static const SolverFactory& factory();

        std::shared_ptr<OptimizationSolver> get_solver(
            const std::string& problem) const;
        inline const std::vector<std::string>& get_solver_names() const
        {
            return solver_names_;
        }

    private:
        SolverFactory();
        std::map<std::string, std::shared_ptr<OptimizationSolver>> solvers_;
        std::vector<std::string> solver_names_;
    };
} // namespace opt
} // namespace ccd
