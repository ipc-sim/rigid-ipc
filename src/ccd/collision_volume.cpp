#include <ccd/collision_volume.hpp>

namespace ccd {

void compute_volumes_fixed_toi(const Eigen::MatrixX2d& V,
    const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
    const EdgeEdgeImpacts& ee_impacts, const Eigen::VectorXi& edge_impact_map,
    const double epsilon, Eigen::VectorXd& volumes)
{

    volumes.resize(E.rows());

    for (long i = 0; i < edge_impact_map.rows(); ++i) {
        if (edge_impact_map[i] == -1) {
            volumes(i) = 0;
            continue;
        }

        EdgeEdgeImpact ee_impact = ee_impacts[size_t(edge_impact_map[i])];

        volumes(i)
            = collision_volume_fixed_toi(V, U, E, ee_impact, int(i), epsilon);
    }
}

double collision_volume_fixed_toi(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpact& impact, const int edge_id, const double epsilon)
{
    double alpha;
    if (edge_id == impact.impacted_edge_index) {
        alpha = impact.impacted_alpha;
    } else if (edge_id == impact.impacting_edge_index) {
        alpha = impact.impacting_alpha;
    } else {
        return 0;
    }
    Eigen::Vector2i e_ij = edges.row(edge_id);

    // we get the position and velocity of the edge vertices
    Eigen::Vector2d Vi = vertices.row(e_ij(0));
    Eigen::Vector2d Vj = vertices.row(e_ij(1));
    Eigen::Vector2d Ui = displacements.row(e_ij(0));
    Eigen::Vector2d Uj = displacements.row(e_ij(1));

    double toi = impact.time;

    return space_time_collision_volume(Vi, Vj, Ui, Uj, toi, alpha, epsilon);
}

double space_time_collision_volume(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const double& toi, const double& alpha,
    const double epsilon)
{

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
    double volume = -(1.0 - toi)
        * std::sqrt(epsilon * epsilon * e_length_toi * e_length_toi
            + U_ij_dot_e_rot90_toi * U_ij_dot_e_rot90_toi);

    return volume;
}

} // namespace ccd
