#include "logger.hpp"

#include <chrono>  // chrono::system_clock
#include <ctime>   // localtime
#include <sstream> // stringstream
#include <iomanip> // put_time
#include <string>  // string

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

    std::string now(){
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y_%m_%d_%X");
        return ss.str();

    }


} // namespace log

} // namespace ccd
