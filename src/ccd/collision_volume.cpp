#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

namespace ccd {

double collision_volume(const Eigen::MatrixX2d& vertices,
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

    return collision_volume(Vi, Vj, Ui, Uj, toi, alpha, epsilon);
}

void collision_volume_grad(const Eigen::MatrixX2d& vertices,
    const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
    Eigen::VectorXd& grad)
{
    ccd::autodiff::ImpactNode impact_node;

    grad.resize(vertices.rows() * 2);
    grad.setZero();

    Eigen::Vector2i e_ij = edges.row(edge_id);
    Eigen::Vector2i e_kl;
    if (edge_id == impact.impacted_edge_index) {
        // impacted edge is this edge, and the impacting node is in
        // the other edge.
        e_kl = edges.row(impact.impacting_edge_index);
        impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vK
                                                   : ccd::autodiff::vL;

    } else if (edge_id == impact.impacting_edge_index) {
        e_kl = edges.row(impact.impacted_edge_index);
        impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vI
                                                   : ccd::autodiff::vJ;
    } else {
        return;
    }

    // we get the position and velocity of the edge vertices
    Eigen::Vector2d Vi = vertices.row(e_ij(0));
    Eigen::Vector2d Vj = vertices.row(e_ij(1));
    Eigen::Vector2d Ui = displacements.row(e_ij(0));
    Eigen::Vector2d Uj = displacements.row(e_ij(1));

    Eigen::Vector2d Vk = vertices.row(e_kl(0));
    Eigen::Vector2d Vl = vertices.row(e_kl(1));
    Eigen::Vector2d Uk = displacements.row(e_kl(0));
    Eigen::Vector2d Ul = displacements.row(e_kl(1));

    Vector8d el_grad = ccd::autodiff::collision_volume_grad(
        Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, impact_node, epsilon);
    // NOTE: local gradient is in format xy, xy, xy, xy
    // while global gradient is in format x,x,x,x ... y, y, y, y
    // TODO: change local gradient format for _consistency_
    auto num_v = vertices.rows();
    grad(e_ij(0)) = el_grad(0);
    grad(num_v + e_ij(0)) = el_grad(1);

    grad(e_ij(1)) = el_grad(2);
    grad(num_v + e_ij(1)) = el_grad(3);

    grad(e_kl(0)) = el_grad(4);
    grad(num_v + e_kl(0)) = el_grad(5);

    grad(e_kl(1)) = el_grad(6);
    grad(num_v + e_kl(1)) = el_grad(7);
}

double collision_volume(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj, const double& toi,
    const double& alpha, const double epsilon)
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

}
