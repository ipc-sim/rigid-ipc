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

#include <ccd/rigid_body/interval.hpp>
#include <utils/eigen_ext.hpp>

namespace ccd {

namespace logger {

    /// @brief Format an eigen MatrixXd
    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 16);

    /// @brief Get a string of the current time
    std::string now();

    /// @brief Format a string for an Interval
    std::string fmt_interval(const Interval& i, const int precision = 16);
    /// @brief Format an eigen VectorX<Interval>
    std::string fmt_eigen_intervals(
        const Eigen::VectorX<Interval>& x, const int precision = 16);

} // namespace logger

} // namespace ccd
