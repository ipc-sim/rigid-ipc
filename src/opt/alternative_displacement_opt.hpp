/**
 * Methods for optimizing the displacments with a penalty constraint.
 */

#include <opt/displacement_opt.hpp>

namespace ccd {
namespace opt {
    namespace alt {

        // Create a OptimizationProblem for displacment optimization
        void setup_displacement_optimization_problem(const Eigen::MatrixX2d& V,
            const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
            const DetectionMethod ccd_detection_method,
            const bool recompute_collision_set, OptimizationProblem& problem);

        // Optimize the displacment opt problem with the given method and
        // starting value.
        OptimizationResults displacement_optimization(
            OptimizationProblem& problem, const Eigen::MatrixX2d& U0,
            SolverSettings& settings);

    } // namespace alt
} // namespace opt
} // namespace ccd
