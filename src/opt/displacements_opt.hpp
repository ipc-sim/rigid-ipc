#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    void displacements_optimization_step(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const double barrier_s,
        const double barrier_beta, const DetectionMethod ccd_detection_method,
        Eigen::MatrixX2d& Uopt);

    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map);

}

}
