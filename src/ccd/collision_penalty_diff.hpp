#pragma once

#include <Eigen/Core>

#include <autodiff/autodiff_types.hpp>
#include <ccd/collision_constraint_diff.hpp>
#include <ccd/collision_detection.hpp>

/**
 * @namespace ccd
 * @brief Continuous collision detection
 */
namespace ccd {

/**
 * @namespace ccd::autodiff
 * @brief CCD functions that use autodiff features
 */
namespace autodiff {

    /**
     * @brief Compute the penalty of an impact between edge_ij and edge_kl.
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

    // ------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Penalties Derivatives
    // ------------------------------------------------------------------------

    /**
     * @brief Compute the value of the collision constraint for all edge-edge
     * impacts.
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] penalties A vector for the computed constraint values.
     */
    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::VectorXd& penalties);

    /**
     * @brief Compute the first derivative of the collision constraint for all
     * edge-edge impacts.
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] constraint_grad  First derivative of the constraint.
     */
    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::MatrixXd& penalty_grad);

    /**
     * @brief Compute the second derivative of the collision constraint for all
     * edge-edge impacts.
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] constraint_grad  Second derivative of the constraint.
     */
    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        std::vector<Eigen::SparseMatrix<double>>& penalty_hessian);

    // ------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Penalties Derivatives
    // ------------------------------------------------------------------------

    /**
     * @brief Compute the value of the collision constraint for a single
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacments        All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              Edge-edge impact
     * @param[in] edge_id             Impact assigned to each edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @return The constraint value for the impact
     */
    double collision_penalty_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon);

    /**
     * @brief Compute the first derivative of the collision constraint for an
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] gradient  First derivative of the constraint.
     */
    void collision_penalty_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::VectorXd& gradient);

    /**
     * @brief Compute the second derivative of the collision constraint for an
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] hessian  Second derivative of the constraint.
     */
    void collision_penalty_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::MatrixXd& hessian);

    /**
     * @brief Compute the first or second derivative of the collision constraint
     * for an edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] derivative  First or second derivative of the constraint.
     */
    template <int I>
    void collision_penalty_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, const int order,
        Eigen::Matrix<double, Eigen::Dynamic, I>& derivative);

    //------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Penalties Derivatives
    //------------------------------------------------------------------------

    /**
     * @brief Compute the derivative of the collision constraint for an
     * edge-edge impact.
     *
     * @param[in] Vi                  First vertex of the impacted edge
     * @param[in] Vj                  Second vertex of the impacted edge
     * @param[in] Vk                  First vertex of the impacting edge
     * @param[in] Vl                  Second vertex of the impacting edge
     * @param[in] Ui                  First displacment of the impacted edge
     * @param[in] Uj                  Second displacment of the impacted edge
     * @param[in] Uk                  First displacment of the impacting edge
     * @param[in] Ul                  Second displacment of the impacting edge
     * @param[in] impact_node         Which node is impacting the other edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @return The differentiable scalar for the constraint
     */
    DScalar collision_penalty_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon);

    /**
     * @brief Compute the derivative of the collision constraint for an
     * edge-edge impact using finite differences.
     *
     * @param[in] Vi                  First vertex of the impacted edge
     * @param[in] Vj                  Second vertex of the impacted edge
     * @param[in] Vk                  First vertex of the impacting edge
     * @param[in] Vl                  Second vertex of the impacting edge
     * @param[in] Ui                  First displacment of the impacted edge
     * @param[in] Uj                  Second displacment of the impacted edge
     * @param[in] Uk                  First displacment of the impacting edge
     * @param[in] Ul                  Second displacment of the impacting edge
     * @param[in] impact_node         Which node is impacting the other edge
     * @param[in] barrier_epsilon     Temporal distance for the start of the
     *                                barrier
     *
     * @param[out] grad  Gradient of the constraint.
     */
    void collision_penalty_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon, Vector8d& grad);

} // namespace autodiff

} // namespace ccd
