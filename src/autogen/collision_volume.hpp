#ifndef CCD_AUTO_COLLISION_VOLUME_HPP
#define CCD_AUTO_COLLISION_VOLUME_HPP

#include <Eigen/Core>
#include <autodiff/autodiff_types.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    using Vector2T = Eigen::Matrix<T, 2, 1>;

    template <typename T>
    void time_of_impact_coeff(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        T& a, T& b, T& c);


    // Computes the space-time volume of intersection between edge_ij and
    // edge_kl
    //
    // Inputs:
    //   V      : vertex positions
    //   U      : vertex displacements
    //
    // Returns:
    //   toi    : time of impact
    template <typename T>
    T _collision_volume_(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        const Vector2T<T>& Ul,
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

    Vector8d collision_volume_grad_fd(
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
