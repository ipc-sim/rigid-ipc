#include <ccd/collision_penalty_diff.hpp>

#include <autodiff/finitediff.hpp>
#include <ccd/time_of_impact.hpp>
#include <opt/barrier.hpp>

#include <algorithm>
#include <iostream>

namespace ccd {
namespace autodiff {

    // Compute the penalty of an impact between edge_ij and edge_kl.
    template <typename T>
    T collision_penalty(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const ImpactNode impact_node,
        const double barrier_epsilon)
    {
        T toi, alpha;
        bool success;

        T t = T(1);                             // Final time for the time step
        T time_scale = t + 2 * barrier_epsilon; // 2x for safety

        const Vector2T<T> scaled_Ui = time_scale * Ui,
                          scaled_Uj = time_scale * Uj,
                          scaled_Uk = time_scale * Uk,
                          scaled_Ul = time_scale * Ul;

        double edge_length = (Vj - Vi).norm();

        switch (impact_node) {
        case vK:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vk, scaled_Ui, scaled_Uj, scaled_Uk, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vk, scaled_Ui, scaled_Uj, scaled_Uk, toi, alpha);
            break;
        case vL:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vl, scaled_Ui, scaled_Uj, scaled_Ul, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vl, scaled_Ui, scaled_Uj, scaled_Ul, toi, alpha);
            break;
        case vI:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vi, scaled_Uk, scaled_Ul, scaled_Ui, toi);
            alpha = T(0);
            break;
        case vJ:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vj, scaled_Uk, scaled_Ul, scaled_Uj, toi);
            alpha = T(1);
            break;
        }

        if (!success) {
            return T(0);
        }

        return edge_length
            * opt::spline_barrier<T>(time_scale * toi - t, barrier_epsilon);
        // return (time_scale * toi - t)
        //     * (scaled_Ui.squaredNorm() + scaled_Uj.squaredNorm()
        //         + scaled_Uk.squaredNorm() + scaled_Ul.squaredNorm());
    }

    // ------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Penalty Derivatives
    // ------------------------------------------------------------------------

    // Compute the value of the collision constraint for all edge-edge impacts.
    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::VectorXd& penalties)
    {
        compute_constraints_refresh_toi(V, U, E, ee_impacts, edge_impact_map,
            barrier_epsilon, collision_penalty<double>, penalties);

        // Use a default value of 10 * barrier_epsilon because any value greater
        // than epsilon results in zero penalty
        // compute_constraints_per_edge_refresh_toi(V, U, E, ee_impacts,
        //     edge_impact_map, barrier_epsilon, collision_penalty<double>,
        //     10 * barrier_epsilon, penalties);
    }

    // Compute the first derivative of the collision constraint for all
    // edge-edge impacts.
    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::MatrixXd& penalties_grad)
    {
        compute_constraints_gradient(V, U, E, ee_impacts, edge_impact_map,
            barrier_epsilon, collision_penalty<DScalar>, penalties_grad);
        // compute_constraints_gradient_per_edge(V, U, E, ee_impacts,
        //     edge_impact_map, barrier_epsilon, collision_penalty<DScalar>,
        //     penalties_grad);
    }

    // Compute the second derivative of the collision constraint for all
    // edge-edge impacts.
    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        std::vector<Eigen::SparseMatrix<double>>& penalties_hessian)
    {
        compute_constraints_hessian(V, U, E, ee_impacts, edge_impact_map,
            barrier_epsilon, collision_penalty<DScalar>, penalties_hessian);
        // compute_constraints_hessian_per_edge(V, U, E, ee_impacts,
        //     edge_impact_map, barrier_epsilon, collision_penalty<DScalar>,
        //     penalties_hessian);
    }

    // ------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Penalties & Derivatives
    // ------------------------------------------------------------------------

    // Compute the value of the collision constraint for a single edge-edge
    // impact.
    double collision_penalty_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon)
    {
        return collision_constraint_refresh_toi(vertices, displacements, edges,
            impact, edge_id, barrier_epsilon, collision_penalty<double>);
    }

    // Compute the first derivative of the collision constraint for an edge-edge
    // impact.
    void collision_penalty_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::VectorXd& gradient)
    {
        Eigen::SparseMatrix<double> sparse_gradient = gradient.sparseView();
        collision_constraint_grad(vertices, displacements, edges, impact,
            edge_id, barrier_epsilon, collision_penalty<DScalar>,
            sparse_gradient);
        gradient = Eigen::VectorXd(sparse_gradient);
    }

    // Compute the second derivative of the collision constraint for an
    // edge-edge impact.
    void collision_penalty_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::MatrixXd& hessian)
    {
        Eigen::SparseMatrix<double> sparse_hessian = hessian.sparseView();
        collision_constraint_hessian(vertices, displacements, edges, impact,
            edge_id, barrier_epsilon, collision_penalty<DScalar>,
            sparse_hessian);
        hessian = Eigen::MatrixXd(sparse_hessian);
    }

    // Compute the first or second derivative of the collision constraint for an
    // edge-edge impact.
    template <int I>
    void collision_penalty_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, const int order,
        Eigen::Matrix<double, Eigen::Dynamic, I>& derivative)
    {
        Eigen::SparseMatrix<double> sparse_derivative = derivative.sparseView();
        collision_constraint_derivative(vertices, displacements, edges, impact,
            edge_id, barrier_epsilon, order, collision_penalty<DScalar>,
            sparse_derivative);
        derivative
            = Eigen::Matrix<double, Eigen::Dynamic, I>(sparse_derivative);
    }

    // ------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Penalty Derivatives
    // ------------------------------------------------------------------------
    // Compute the derivative of the collision constraint for an edge-edge
    // impact.
    DScalar collision_penalty_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon)
    {
        return collision_constraint_differentiable(Vi, Vj, Vk, Vl, Ui, Uj, Uk,
            Ul, impact_node, barrier_epsilon, collision_penalty<DScalar>);
    }

    // Compute the derivative of the collision constraint for an edge-edge
    // impact using finite differences.
    void collision_penalty_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon, Vector8d& grad)
    {
        collision_constraint_grad_fd(Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul,
            impact_node, barrier_epsilon, collision_penalty<double>, grad);
    }

    template double collision_penalty<double>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, const ImpactNode, const double);

    template DScalar collision_penalty<DScalar>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        DVector2 const&, DVector2 const&, DVector2 const&, DVector2 const&,
        const ImpactNode, const double);

    template void collision_penalty_derivative<1>(const Eigen::MatrixX2d&,
        const Eigen::MatrixX2d&, const Eigen::MatrixX2i&, const EdgeEdgeImpact&,
        const int, const double, const int, Eigen::VectorXd&);

    template void collision_penalty_derivative<Eigen::Dynamic>(
        const Eigen::MatrixX2d&, const Eigen::MatrixX2d&,
        const Eigen::MatrixX2i&, const EdgeEdgeImpact&, const int, const double,
        const int, Eigen::MatrixXd&);

} // namespace autodiff
} // namespace ccd
