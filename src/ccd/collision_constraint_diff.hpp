/**
 * General code to compute the collision constraint and its derivative(s).
 */
#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <autodiff/autodiff_types.hpp>
#include <ccd/collision_detection.hpp>

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

    /// @brief Generic type 2×1 vector
    template <typename T> using Vector2T = Eigen::Matrix<T, 2, 1>;

    /// @brief Enumeration of the for end points of a edge-edge impact
    enum ImpactNode {
        vI, ///< @brief First vertex of impacted edge
        vJ, ///< @brief Second vertex of impacted edge
        vK, ///< @brief First vertex of impacting edge
        vL  ///< @brief Second vertex of impacting edge
    };

    /// @brief Function for computing the constraints of a single impact
    template <typename T>
    using constraint_func = std::function<T(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Vector2T<T>& Ui, const Vector2T<T>& Uj,
        const Vector2T<T>& Uk, const Vector2T<T>& Ul,
        const ImpactNode impact_node, const double epsilon)>;

    ///////////////////////////////////////////////////////////////////////////
    // Global constraints over the entire simulation

    // One constraint per edge versions:

    /**
     * @brief Compute the constraints globally where each edge has a single
     * constraint.
     *
     * @param[in] V                   All vertices positions.
     * @param[in] U                   All vertices displacements.
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges.
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     * @param[in] default_value       The value used when there is no collision
     *
     * @param[out] constraints A vector for the computed constraint values
     */
    void compute_constraints_per_edge_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<double>& compute_constraint,
        const double default_value, Eigen::VectorXd& constraints);

    /**
     * @brief Compute the first derivative of the constraints globally where
     * each edge has a single constraint.
     *
     * @param[in] V                   All vertices positions.
     * @param[in] U                   All vertices displacements.
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges.
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] constraints_grad A matrix for the computed gradient of
     *                              constraint values
     */
    void compute_constraints_gradient_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraints_grad);

    /**
     * @brief Compute the second derivative of the constraints globally where
     * each edge has a single constraint.
     *
     * @param[in] V                   All vertices positions.
     * @param[in] U                   All vertices displacements.
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges.
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] constraints A std::vector of matrices for the computed
     *                         hessian of the constraint values
     */
    void compute_constraints_hessian_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        std::vector<Eigen::MatrixXd>& constraints_hessian);

    // One constraint per impact versions:

    /**
     * @brief Compute the value of the collision constraint for all edge-edge
     * impacts
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] constraints A vector for the computed constraint values.
     */
    void compute_constraints_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints);

    void compute_constraints_dense_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints);

    /**
     * @brief Compute the first derivative of the collision constraint for all
     * edge-edge impacts.
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] constraint_grad  First derivative of the constraint.
     */
    void compute_constraints_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraint_grad);

    void compute_constraints_dense_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraint_grad);

    /**
     * @breif Compute the second derivative of the collision constraint for all
     * edge-edge impacts.
     *
     * @param[in] V                   All vertices positions
     * @param[in] U                   All vertices displacements
     * @param[in] E                   Edges as pair of vertex indices
     * @param[in] ee_impacts          List of impact between two edges
     * @param[in] edge_impact_map     Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] constraint_grad  Vector of Matricies of the second derivative
     *                              of the constraint
     */
    void compute_constraints_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        std::vector<Eigen::SparseMatrix<double>>& constraint_hessian);

    ///////////////////////////////////////////////////////////////////////////
    // SINGLE IMPACT GLOBAL Constraints Derivatives

    /**
     * @brief Compute the value of the collision constraint for a single
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacments        All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              Edge-edge impact
     * @param[in] edge_id             Impact assigned to each edge
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @return The constraint value for the impact
     */
    double collision_constraint_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<double>& compute_constraint);

    /**
     * @brief Compute the first derivative of the collision constraint for an
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] gradient  First derivative of the constraint.
     */
    void collision_constraint_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& gradient);

    /**
     * @brief Compute the second derivative of the collision constraint for an
     * edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] hessian  Second derivative of the constraint.
     */
    void collision_constraint_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& hessian);

    /**
     * @brief Compute the first or second derivative of the collision constraint
     * for an edge-edge impact.
     *
     * @param[in] vertices            All vertices positions.
     * @param[in] displacements       All vertices displacements.
     * @param[in] edges               Edges as pair of vertex indices
     * @param[in] impact              An impact between two edges.
     * @param[in] edge_id             The edge for which constraint is computed
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] derivative  First or second derivative of the constraint.
     */
    void collision_constraint_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const int order, const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& derivative);

    ////////////////////////////////////////////////////////////////////////////
    // Single impact local Constraints Derivatives

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
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @return The differentiable scalar for the constraint
     */
    DScalar collision_constraint_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon,
        const constraint_func<DScalar>& compute_constraint);

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
     * @param[in] epsilon             A number used during computation
     *                                (e.g. volume epsilon or barrier epsilon)
     * @param[in] compute_constraint  A function to compute the constraint
     *                                locally
     *
     * @param[out] grad  Gradient of the constraint.
     */
    void collision_constraint_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon, const constraint_func<double>& compute_constraint,
        Vector8d& grad);

} // namespace autodiff

} // namespace ccd
