#pragma once

#include <Eigen/Core>

#include <autodiff/autodiff_types.hpp>
#include <ccd/collision_detection.hpp>
#include <ccd/collision_volume_diff.hpp>

/**
 * @namespace ccd:
 * @brief
 */
namespace ccd {

/**
 * @namespace ccd::autodiff
 * @brief CCD functions that use autodiff features
 */
namespace autodiff {
    // ------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Penalties Derivatives
    // ------------------------------------------------------------------------
    /**
     * Compute the first derivative of the collision penalty for all edge-edge
     * impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] ee_impacts       : List of impact between two edges.
     *  @param[in] edge_impact_map  : Impact assigned to each edge
     *
     *  @param[out] derivative      : Matrix of first derivative of the penalty
     *                                size (2*V, E)
     *
     */
    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::MatrixXd& penalty_grad);

    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        std::vector<Eigen::MatrixXd>& penalty_hessian);

    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::VectorXd& penalties);

    // ------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Penalties Derivatives
    // ------------------------------------------------------------------------
    double collision_penalty_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon);
    /**
     * Compute the first or second derivative of the collision penalty for an
     * edge-edge impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] edge_id          : The edge for which penalty is computed
     *
     *  @param[out] derivative      : first or second derivative of the
     * penalty.
     */
    template <int T>
    void collision_penalty_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, const int order,
        Eigen::Matrix<double, Eigen::Dynamic, T>& derivative);

    void collision_penalty_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::VectorXd& gradient);

    void collision_penalty_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::MatrixXd& hessian);

    // ------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Penalties Derivatives
    // ------------------------------------------------------------------------

    /**
     * Compute the penalty of intersection between for edge_ij for an impact
     * with edge_kl
     *
     *  @param V_{ijkl}     : Vertices positions.
     *  @param U_{ijkl}     : Vertices displacements.
     *  @param impact_node  : The node (i,j,k or l) that caused the impact
     *
     *  @return             : The space-time interference penalty.
     */
    template <typename T>
    T collision_penalty(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const ImpactNode impact_node,
        const double barrier_epsilon);

    DScalar collision_penalty_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon);

    ///@brief helper function for testing
    void collision_penalty_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon, Vector8d& grad);

} // namespace autodiff

} // namespace ccd
