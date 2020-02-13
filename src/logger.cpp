#include "logger.hpp"

#include <chrono>  // chrono::system_clock
#include <ctime>   // localtime
#include <iomanip> // put_time
#include <sstream> // stringstream
#include <string>  // string

#include <spdlog/sinks/stdout_color_sinks.h>

namespace ccd {
namespace logger {

    static const Eigen::IOFormat
        CleanFmt(Eigen::FullPrecision, 0, ", ", "\n", "", "");
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

    std::string fmt_interval(const Interval& i, const int precision)
    {
        std::stringstream ssx;
        ssx.precision(precision);
        ssx << "[" << i.lower() << ", " << i.upper() << "]";
        return ssx.str();
    }

    std::string
    fmt_eigen_intervals(const Eigen::VectorX<Interval>& x, const int precision)
    {
        std::stringstream ss;
        ss << "[";
        for (int i = 0; i < x.size(); i++) {
            ss << fmt_interval(x(i)) << (i < x.size() - 1 ? ", " : "");
        }
        ss << "]";
        return ss.str();
    }

} // namespace logger
} // namespace ccd
