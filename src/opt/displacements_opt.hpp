#pragma once

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    enum OptimizationMethod {
        MMA,   ///<@brief Method of Moving Asymptotes (NLopt)
        SLSQP, ///<@brief Sequential Least-Squares Quadratic Programming (NLopt)
        IP     ///<@brief Interior-Point Method (Ipopt)
    };

    void displacements_optimization(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double volume_epsilon, const DetectionMethod ccd_detection_method,
        const OptimizationMethod opt_method, const unsigned max_iter,
        Eigen::MatrixX2d& Uopt);

    void detect_collisions(const Eigen::MatrixX2d& V, const Eigen::MatrixX2d& U,
        const Eigen::MatrixX2i& E, const ccd::DetectionMethod detection_method,
        EdgeEdgeImpacts& ee_impacts, Eigen::VectorXi& edge_impact_map);

}

}
