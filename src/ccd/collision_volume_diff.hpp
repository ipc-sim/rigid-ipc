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

    /**
    * Compute the volume of intersection between a edge and a vertex
    *
    *
    *     @param V_{ijk}          : Vertices positions.
    *     @param U_{ijk}          : Vertices displacements.
    *     @param epsilon          : The time scale used for minimal volume.
    *
    * @return                     : The space-time interference volume.
    */
    template <typename T>
    T collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        const double epsilon);


    Vector8d collision_volume_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const double epsilon);

    Vector8d collision_volume_grad_fd(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const double epsilon);
}

}
#endif // COLLISION_VOLUME_DIFF_H
