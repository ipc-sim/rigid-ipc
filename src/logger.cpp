#include "logger.hpp"

#include <ctime>   // localtime
#include <iomanip> // put_time
#include <sstream> // stringstream
#include <string>  // string

#include <spdlog/fmt/chrono.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <ipc/utils/logger.hpp>

#include <tbb/parallel_sort.h>

#include <interval/interval.hpp>

namespace ipc::rigid {

static const Eigen::IOFormat RowVectorFmt(
    Eigen::FullPrecision, Eigen::DontAlignCols, ", ", ",\n", "[", "]", "", "");
static const Eigen::IOFormat MatrixFmt(
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
    ssx << std::setprecision(precision);
    if (m.rows() == 1) {
        ssx << m.format(RowVectorFmt);
    } else {
        ssx << m.format(MatrixFmt);
    }
    return ssx.str();
}

std::string current_time_string()
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

std::string fmt_eigen_intervals(const VectorXI& x, const int precision)
{
    std::stringstream ss;
    ss << std::setprecision(precision) << "[";
    for (int i = 0; i < x.size(); i++) {
        ss << fmt_interval(x(i)) << (i < x.size() - 1 ? ", " : "");
    }
    ss << "]";
    return ss.str();
}

void set_logger_level(spdlog::level::level_enum log_level)
{
    spdlog::set_level(log_level);
    ipc::logger().set_level(log_level);
    // fd::logge().set_logger_level(log_level);
}

void print_candidates(const Candidates& candidates)
{
    fmt::print("fv_candidates: [");
    for (const auto& fv_candidate : candidates.fv_candidates) {
        fmt::print(
            "({:d}, {:d}), ", fv_candidate.face_index,
            fv_candidate.vertex_index);
    }
    fmt::print("]\nee_candidates: [");
    for (const auto& ee_candidate : candidates.ee_candidates) {
        fmt::print(
            "({:d}, {:d}), ", ee_candidate.edge0_index,
            ee_candidate.edge1_index);
    }
    fmt::print("]\n");
}

void print_candidates_sorted(Candidates candidates)
{
    tbb::parallel_sort(
        candidates.ev_candidates.begin(), candidates.ev_candidates.end(),
        [](const EdgeVertexCandidate& ev0, const EdgeVertexCandidate& ev1) {
            if (ev0.edge_index == ev1.edge_index) {
                return ev0.vertex_index < ev1.vertex_index;
            }
            return ev0.edge_index < ev1.edge_index;
        });

    tbb::parallel_sort(
        candidates.ee_candidates.begin(), candidates.ee_candidates.end(),
        [](const EdgeEdgeCandidate& ee0, const EdgeEdgeCandidate& ee1) {
            size_t e0_min = std::min(ee0.edge0_index, ee0.edge1_index);
            size_t e0_max = std::max(ee0.edge0_index, ee0.edge1_index);
            size_t e1_min = std::min(ee1.edge0_index, ee1.edge1_index);
            size_t e1_max = std::max(ee1.edge0_index, ee1.edge1_index);
            if (e0_min == e1_min) {
                return e0_max < e1_max;
            }
            return e0_min < e1_min;
        });

    tbb::parallel_sort(
        candidates.fv_candidates.begin(), candidates.fv_candidates.end(),
        [](const FaceVertexCandidate& fv0, const FaceVertexCandidate& fv1) {
            if (fv0.face_index == fv1.face_index) {
                return fv0.vertex_index < fv1.vertex_index;
            }
            return fv0.face_index < fv1.face_index;
        });

    print_candidates(candidates);
}

} // namespace ipc::rigid
