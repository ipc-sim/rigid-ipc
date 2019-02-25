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
        // IJ is the impactED edge. The impacting node is K or L.
        e_kl = edges.row(impact.impacting_edge_index);
        impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vK
                                                   : ccd::autodiff::vL;

    } else if (edge_id == impact.impacting_edge_index) {
        // IJ is the impactING edge
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

    // LOCAL gradient and hessian. Indices refer to the 4 vertices involded
    // in the collision
    Vector8d el_grad;
    Matrix8d el_hessian;

    ccd::autodiff::collision_volume_grad(
        Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, impact_node, epsilon, el_grad, el_hessian);

    // Assemble into GLOBAL gradient and hessian
    auto num_vertices = vertices.rows();
    int nodes[4] = { e_ij(0), e_ij(1), e_kl(0), e_kl(1) };
    // Note: global gradient is sorted as x,x,x,...y,y,y
    // while local gradient is sorted as x,y,x,y,...,x,y
    for (int i = 0; i < 4; i++) {
        grad(nodes[i]) = el_grad(2 * i);
        grad(nodes[i] + num_vertices) = el_grad(2 * i + 1);
    }
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
