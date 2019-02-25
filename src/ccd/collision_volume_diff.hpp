#ifndef COLLISION_VOLUME_DIFF_H
#define COLLISION_VOLUME_DIFF_H

#include <Eigen/Core>

#include <autodiff/autodiff_types.hpp>

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
    template <typename T>
    using Vector2T = Eigen::Matrix<T, 2, 1>;

    enum ImpactNode {
        vI,
        vJ,
        vK,
        vL
    };

    /**
    * Compute the volume of intersection between for edge_ij for an impact with
    * edge_kl
    *
    *
    *     @param V_{ijkl}          : Vertices positions.
    *     @param U_{ijkl}          : Vertices displacements.
    *     @param impact_node       : The node (i,j,k or l) that caused the impact
    *     @param epsilon           : The time scale used for minimal volume.
    *
    * @return                     : The space-time interference volume.
    */
    template <typename T>
    T collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        const Vector2T<T>& Ul,
        const ImpactNode impact_node,
        const double epsilon);

    void collision_volume_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const ImpactNode impact_node,
        const double epsilon,
        Vector8d& grad,
        Matrix8d& hessian);

    void collision_volume_grad_fd(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const ImpactNode impact_node,
        const double epsilon,
        Vector8d& grad);

}

}
#endif // COLLISION_VOLUME_DIFF_H
