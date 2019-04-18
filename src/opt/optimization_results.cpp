#include <opt/optimization_results.hpp>

namespace ccd {
namespace opt {

    OptimizationResults::OptimizationResults()
        : x()
        , minf(std::numeric_limits<double>::infinity())
        , success()
    {
    }

    OptimizationResults::OptimizationResults(
        Eigen::MatrixXd x, double minf, bool success)
        : x(x)
        , minf(minf)
        , success(success)
    {
    }

} // namespace opt
} // namespace ccd
