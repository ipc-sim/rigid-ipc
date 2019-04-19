#pragma once

#include <Eigen/Core>

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
    template <typename T> using Vector2T = Eigen::Matrix<T, 2, 1>;

    enum ImpactNode { vI, vJ, vK, vL };

    template <typename T>
    using constraint_func = std::function<T(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Vector2T<T>& Ui, const Vector2T<T>& Uj,
        const Vector2T<T>& Uk, const Vector2T<T>& Ul,
        const ImpactNode impact_node, const double epsilon)>;

    // -----------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Constraints Derivatives
    // -----------------------------------------------------------------------------

    void compute_constraints_per_edge_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints);
    void compute_constraints_gradient_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraints_grad);
    void compute_constraints_hessian_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        std::vector<Eigen::MatrixXd>& constraints_hessian);

    void compute_constraints_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints);

    /**
     * Compute the first derivative of the collision constraint for all
     * edge-edge impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] ee_impacts       : List of impact between two edges.
     *  @param[in] edge_impact_map  : Impact assigned to each edge
     *  @param[in] epsilon          : The time scale used for minimal
     *                                constraint.
     *
     *  @param[out] derivative      : Matrix of first derivative of the
     *                                constraint size (2*V, E)
     */
    void compute_constraints_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,

        Eigen::MatrixXd& constraint_grad);

    void compute_constraints_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,

        std::vector<Eigen::MatrixXd>& constraint_hessian);

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Constraints Derivatives
    // -----------------------------------------------------------------------------
    double collision_constraint_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<double>& compute_constraint);
    /**
     * Compute the first or second derivative of the collision constraint for an
     * edge-edge impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] edge_id          : The edge for which constraint is computed
     *  @param[in] epsilon          : The time scale used for minimal
     *                                constraint.
     *
     *  @param[out] derivative      : first or second derivative of the
     *                                constraint.
     */
    template <int I>
    void collision_constraint_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const int order, const constraint_func<DScalar>& compute_constraint,

        Eigen::Matrix<double, Eigen::Dynamic, I>& derivative);

    void collision_constraint_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::VectorXd& gradient);

    void collision_constraint_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& hessian);

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Constraints Derivatives
    // -----------------------------------------------------------------------------

    DScalar collision_constraint_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon,
        const constraint_func<DScalar>& compute_constraint);

    ///@brief helper function for testing
    void collision_constraint_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon, const constraint_func<double>& compute_constraint,
        Vector8d& grad);

} // namespace autodiff

} // namespace ccd
