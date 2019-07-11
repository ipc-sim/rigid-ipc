#include "optimization_solver.hpp"

namespace ccd {
namespace opt {

    OptimizationSolver::OptimizationSolver()
        : free_dof()
        , max_iterations(3000)
    {
    }
    OptimizationSolver::OptimizationSolver(const int max_iterations)
        : free_dof()
        , max_iterations(max_iterations)
    {
    }

    OptimizationSolver::~OptimizationSolver() {}

    // Initialize free_dof with indices of dof that are not fixed.
    void OptimizationSolver::init_free_dof(Eigen::VectorXb is_dof_fixed)
    {
        free_dof = Eigen::VectorXi(is_dof_fixed.size() - is_dof_fixed.count());
        for (int i = 0, j = 0; i < is_dof_fixed.size(); i++) {
            if (!is_dof_fixed(i)) {
                free_dof(j++) = i;
            }
        }
    }

} // namespace opt
} // namespace ccd
