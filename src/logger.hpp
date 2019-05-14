#pragma once

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#endif

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

#include <Eigen/Core>

namespace ccd {

namespace log {

    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 10);
    std::string now();

} // namespace log

} // namespace ccd
