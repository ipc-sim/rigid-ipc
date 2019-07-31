#include "logger.hpp"

#include <chrono>  // chrono::system_clock
#include <ctime>   // localtime
#include <iomanip> // put_time
#include <sstream> // stringstream
#include <string>  // string

#include <spdlog/sinks/stdout_color_sinks.h>

namespace ccd {
namespace log {
    static const Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "", "");
    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision)
    {
        std::stringstream ssx;
        Eigen::MatrixXd m = x;
        if (m.cols() == 1) {
            m.transposeInPlace();
        }
        ssx << std::setprecision(precision) << m.format(CleanFmt);
        return ssx.str();
    }

    std::string now()
    {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);

        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y_%m_%d_%EX");
        return ss.str();
    }



} // namespace log

} // namespace ccd
