#ifndef CCD_COLLISION_VOLUME_H
#define CCD_COLLISION_VOLUME_H

#include <Eigen/Core>
#include <FixingCollisions/collision_detection.hpp>

namespace ccd {
double collision_volume(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges,
    const Impact& impact,
    const Eigen::VectorXd grad_volume);
}

#endif
