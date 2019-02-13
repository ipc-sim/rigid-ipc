#ifndef CCD_COLLISION_VOLUME_H
#define CCD_COLLISION_VOLUME_H

#include <Eigen/Core>
#include <FixingCollisions/collision_detection.hpp>

namespace ccd {

/**
Compute the volume of the STIV.

\f$V = (1-\tau_I)\sqrt{\epsilon^2 \|e(\tau_I)\|^2 + (U_{ij} \cdot
e(\tau_I)^\perp)^2}\f$

@param vertices Matrix of vertices positions.
@param displacements Matrix of displacment vectors.
@param edges Matrix of edges as vertex indices.
@param impact An impact between two edges.
@param epsilon The time scale used for minimal volume.
@param grad_volue Matrix to store the gradient of the volume in.
@return The space-time interference volume.
*/
double collision_volume(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpactPtr impact, const double epsilon,
    Eigen::MatrixX2d& grad_volume);

}

#endif
