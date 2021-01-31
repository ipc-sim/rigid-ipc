#pragma once

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <utils/eigen_ext.hpp>

#include <ipc/spatial_hash/collision_candidate.hpp>

namespace ipc::rigid {

/// @brief Format an eigen MatrixXd
std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 16);

/// @brief Get a string of the current time
std::string current_time_string();

// NOTE: Interval logging is defined in the CPP file but declares in
// interval.hpp.

void set_logger_level(spdlog::level::level_enum log_level);

void print_candidates(const Candidates& candidates);
void print_candidates_sorted(Candidates candidates);

} // namespace ipc::rigid
