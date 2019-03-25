#include "logger.hpp"

namespace ccd {
namespace log {
    static const Eigen::IOFormat CommaFmt(Eigen::StreamPrecision,
        Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision)
    {
        std::stringstream ssx;
        ssx << std::setprecision(precision) << x.format(CommaFmt);
        return ssx.str();
    }

} // namespace log

} // namespace ccd
