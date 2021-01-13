#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <vector>

namespace swr {

bool write_png(
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& R,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& G,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& B,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& A,
    const std::string& filename);

} // namespace swr
