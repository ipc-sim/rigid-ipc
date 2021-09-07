#pragma once

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>
#include <vector>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <utils/eigen_ext.hpp>

#include <ipc/broad_phase/collision_candidate.hpp>

namespace ipc::rigid {

/// @brief Format an eigen MatrixXd
std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 16);

/// @brief Get a string of the current time
std::string current_time_string();

// NOTE: Interval logging is defined in the CPP file but declares in
// interval.hpp.

const std::vector<std::pair<std::string, spdlog::level::level_enum>>
    SPDLOG_LEVEL_NAMES_TO_LEVELS = { { "trace", spdlog::level::trace },
                                     { "debug", spdlog::level::debug },
                                     { "info", spdlog::level::info },
                                     { "warning", spdlog::level::warn },
                                     { "error", spdlog::level::err },
                                     { "critical", spdlog::level::critical },
                                     { "off", spdlog::level::off } };

void set_logger_level(spdlog::level::level_enum log_level);

void print_candidates(const Candidates& candidates);
void print_candidates_sorted(Candidates candidates);

} // namespace ipc::rigid
