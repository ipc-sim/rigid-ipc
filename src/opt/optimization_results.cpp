#include <opt/optimization_results.hpp>

namespace ccd {
namespace opt {

    OptimizationResults::OptimizationResults()
        : x()
        , minf(std::numeric_limits<double>::infinity())
        , success(false)
        , finished(false)
    {
    }

    OptimizationResults::OptimizationResults(
        Eigen::MatrixXd x, double minf, bool success, bool finished)
        : x(x)
        , minf(minf)
        , success(success)
        , finished(finished)
    {
    }

} // namespace opt
} // namespace ccd
