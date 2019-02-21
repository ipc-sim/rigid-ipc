#ifndef CCD_COLLISION_VOLUME_H
#define CCD_COLLISION_VOLUME_H

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>

namespace ccd {

/**
* Compute the collision volume for an edge-edge impact
*
*     @param vertices         : All vertices positions.
*     @param displacements    : All vertices displacements.
*     @param edges            : Edges as pair of vertex indices
*     @param impact           : An impact between two edges.
*     @param edge_id          : The edge for which we are computing the volume
*     @param epsilon          : The time scale used for minimal volume.
*
* @return                     : The space-time interference volume.
*/
double collision_volume(
    const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements,
    const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpactPtr impact,
    const int edge_id,
    const double epsilon);

/**
* Compute the volume of intersection for an edge given a time
* of intersection (toi) and position of intersection (alpha)
*
* \f$V = (1-\tau_I)\sqrt{\epsilon^2 \|e(\tau_I)\|^2 + (U_{ij} \cdot
* e(\tau_I)^\perp)^2}\f$
*
*     @param V_{ij}          : Vertices positions.
*     @param U_{ij}          : Vertices displacements.
*     @param toi             : Time of impact.
*     @param alpha           : Position along edge where impact happens.
*
* @return                     : The space-time interference volume.
*/
double collision_volume(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const double& toi,
    const double& alpha,
    const double epsilon);
}

#endif
