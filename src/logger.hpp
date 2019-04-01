#pragma once

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <Eigen/Core>

namespace ccd {

namespace log {

    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 10);
    std::string now();

} // namespace log

} // namespace ccd
