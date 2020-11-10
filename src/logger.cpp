#include "logger.hpp"

#include <ctime> // localtime
#include <fmt/chrono.h>
#include <iomanip> // put_time
#include <sstream> // stringstream
#include <string>  // string

#include <spdlog/sinks/stdout_color_sinks.h>

#include <ipc/utils/logger.hpp>

#include <interval/interval.hpp>

namespace ccd {
namespace logger {

    static const Eigen::IOFormat CleanFmt(
        Eigen::FullPrecision,
        Eigen::DontAlignCols,
        ", ",
        ",\n",
        "[",
        "]",
        "[",
        "]");
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
        return fmt::format("{:%F-%H-%M-%S}", fmt::localtime(std::time(0)));
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
        ss << std::setprecision(precision) << "[";
        for (int i = 0; i < x.size(); i++) {
            ss << fmt_interval(x(i)) << (i < x.size() - 1 ? ", " : "");
        }
        ss << "]";
        return ss.str();
    }

    void set_level(spdlog::level::level_enum log_level)
    {
        spdlog::set_level(log_level);
        ipc::logger().set_level(log_level);
        // fd::logger().set_level(log_level);
    }

} // namespace logger
} // namespace ccd
