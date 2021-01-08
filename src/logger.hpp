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

#include <utils/eigen_ext.hpp>

#include <ipc/spatial_hash/collision_candidate.hpp>

namespace ccd {

namespace logger {

    /// @brief Format an eigen MatrixXd
    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 16);

    /// @brief Get a string of the current time
    std::string now();

    // NOTE: Interval logging is defined in the CPP file but declares in
    // interval.hpp.

    void set_level(spdlog::level::level_enum log_level);

    void print_candidates(const ipc::Candidates& candidates);
    void print_candidates_sorted(ipc::Candidates candidates);

} // namespace logger

} // namespace ccd
