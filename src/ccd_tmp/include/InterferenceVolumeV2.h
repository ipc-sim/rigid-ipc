#ifndef IAGM_INTERFERENCE_VOLUME_V2_H
#define IAGM_INTERFERENCE_VOLUME_V2_H

#include <Eigen/Core>

namespace IAGM {

// Computes the intersection volume for a edge-vertex pair.
// x0: position of vertex and edge nodes at time t=0
// x1: position of vertex and edge nodes at time t=1
double interference_volume(
    const Eigen::Vector2d x0[3],
    const Eigen::Vector2d x1[3],
    const double min_sep, const double tolerance,
    double toi = -1.0, double alpha = -1.0);
}

bool time_of_collision(
        const Eigen::Vector2d x00, const Eigen::Vector2d x01,
        const Eigen::Vector2d x10, const Eigen::Vector2d x11,
        const Eigen::Vector2d x20, const Eigen::Vector2d x21,
        const double min_sep, const double tolerance,
        double &toi, double &alpha);
#endif
