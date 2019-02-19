#ifndef CCD_AUTO_TOI_HPP
#define CCD_AUTO_TOI_HPP

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

}
}
#endif
