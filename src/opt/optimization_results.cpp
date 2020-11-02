#include <opt/optimization_results.hpp>

namespace ccd {
namespace opt {

    OptimizationResults::OptimizationResults()
        : x()
        , minf(std::numeric_limits<double>::infinity())
        , success(false)
        , finished(false)
        , num_iterations(0)
    {
    }

    OptimizationResults::OptimizationResults(
        Eigen::MatrixXd x,
        double minf,
        bool success,
        bool finished,
        int num_iterations)
        : x(x)
        , minf(minf)
        , success(success)
        , finished(finished)
        , num_iterations(num_iterations)
    {
    }

} // namespace opt
} // namespace ccd
