#ifndef CCD_AUTO_COLLISION_VOLUME_HPP
#define CCD_AUTO_COLLISION_VOLUME_HPP

#include <Eigen/Core>
#include <autodiff/autodiff_types.hpp>

namespace ccd {
namespace autogen {

    template<typename T>
    T _collision_volume_(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Matrix<T, 2, 1>& Ui,
        const Eigen::Matrix<T, 2, 1>& Uj,
        const Eigen::Matrix<T, 2, 1>& Uk,
        const Eigen::Matrix<T, 2, 1>& Ul,
        const double epsilon);

    double collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon);

    Vector8d collision_volume_grad(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul,
        const double epsilon);

}
}
#endif
