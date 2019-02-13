#include <FixingCollisions/collision_volume.hpp>

namespace ccd {

double collision_volume(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpactPtr impact, const double epsilon,
    Eigen::MatrixX2d& /*grad_volume*/)
{
    Eigen::Vector2i e_ij = edges.row(impact->impacted_edge_index);

    // we get the position and velocity of the edge vertices
    Eigen::Vector2d Vi = vertices.row(e_ij(0));
    Eigen::Vector2d Vj = vertices.row(e_ij(1));
    Eigen::Vector2d Ui = displacements.row(e_ij(0));
    Eigen::Vector2d Uj = displacements.row(e_ij(1));

    double toi = impact->time;
    double alpha = impact->impacted_alpha;
    assert(toi >= 0.0 && toi <= 1.0);
    assert(alpha >= 0.0 && alpha <= 1.0);

    // get edge at time of impact (toi)
    Eigen::Vector2d e_toi = (Vj + toi * Uj) - (Vi + toi * Ui);
    Eigen::Vector2d e_rot90_toi(e_toi(1), -e_toi(0));
    double e_length_toi = e_toi.norm();
    assert(e_length_toi > 0.0);

    // get velocity of point of contact along the edge
    Eigen::Vector2d U_ij = Ui + alpha * (Uj - Ui);

    // if velocity is perp. to normal, we must use a non-zero epsilon
    double U_ij_dot_e_rot90_toi = U_ij.dot(e_rot90_toi);
    assert(epsilon > 0 || std::abs(U_ij_dot_e_rot90_toi) > 0);

    // volume = (1-t)\sqrt{\epsilon^2 \|e(t)\|^2 + (U_{ij} \cdot e(t)^\perp)^2}
    double volume = (1.0 - toi)
        * std::sqrt(epsilon * epsilon * e_length_toi * e_length_toi
            + U_ij_dot_e_rot90_toi * U_ij_dot_e_rot90_toi);

    // Flip the volume if we used the wrong perpendicular vector to e_toi
    volume *= volume > 0 ? -1 : 1;

    // TODO: compute gradient
    return volume;
}
}
