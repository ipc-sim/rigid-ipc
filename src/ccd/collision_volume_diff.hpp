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

    // -----------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Volumes Derivatives
    // -----------------------------------------------------------------------------
    /**
     * Compute the first derivative of the collision volume for all edge-edge
     * impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] ee_impacts       : List of impact between two edges.
     *  @param[in] edge_impact_map  : Impact assigned to each edge
     *  @param[in] epsilon          : The time scale used for minimal volume.
     *
     *  @param[out] derivative      : Matrix of first derivative of the volume
     *                                size (2*V, E)
     *
     */
    void compute_volumes_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::MatrixXd& volume_grad);

    void compute_volumes_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        std::vector<Eigen::MatrixXd>& volume_hessian);

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Volumes Derivatives
    // -----------------------------------------------------------------------------

    /**
     * Compute the first or second derivative of the collision volume for an
     * edge-edge impact
     *
     *  @param[in] vertices         : All vertices positions.
     *  @param[in] displacements    : All vertices displacements.
     *  @param[in] edges            : Edges as pair of vertex indices
     *  @param[in] impact           : An impact between two edges.
     *  @param[in] edge_id          : The edge for which volume is computed
     *  @param[in] epsilon          : The time scale used for minimal volume.
     *
     *  @param[out] derivative      : first or second derivative of the
     * volume.
     */
    template <int T>
    void collision_volume_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const int order, Eigen::Matrix<double, Eigen::Dynamic, T>& derivative);

    void collision_volume_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        Eigen::VectorXd& gradient);

    void collision_volume_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        Eigen::MatrixXd& hessian);

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Volumes Derivatives
    // -----------------------------------------------------------------------------

    /**
     * Compute the volume of intersection between for edge_ij for an impact with
     * edge_kl
     *
     *  @param V_{ijkl}     : Vertices positions.
     *  @param U_{ijkl}     : Vertices displacements.
     *  @param impact_node  : The node (i,j,k or l) that caused the impact
     *  @param epsilon      : The time scale used for minimal volume.
     *
     *  @return             : The space-time interference volume.
     */
    template <typename T>
    T collision_volume(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const ImpactNode impact_node,
        const double epsilon);

    DScalar collision_volume_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon);

    ///@brief helper function for testing
    void collision_volume_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon, Vector8d& grad);

}

}
