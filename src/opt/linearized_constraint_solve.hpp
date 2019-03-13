#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    void linearized_constraint_solve(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const DetectionMethod CCD_DETECTION_METHOD,
        const double VOLUME_EPSILON);

}

} // namespace ccd
