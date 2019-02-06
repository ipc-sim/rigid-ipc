#ifndef CCD_COLLISION_VOLUME_H
#define CCD_COLLISION_VOLUME_H

#include <Eigen/Core>
#include <FixingCollisions/collision_detection.hpp>

namespace ccd {
double collision_volume(
    const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const Eigen::MatrixX2i& edges,
    const Impact& impact,
    const Eigen::VectorXd grad_volume);
}

#endif
