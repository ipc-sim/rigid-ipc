#include "time_of_impact_coeff.hpp"

#include <iostream>

#include <autodiff/autodiff_types.hpp>
#include <autodiff/finitediff.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    void time_of_impact_coeff(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        T& a, T& b, T& c)
    {
        // {{toi_abc_ccode}}
    }

    template void time_of_impact_coeff<double>(
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        double&, double&, double&);

    template void time_of_impact_coeff<DScalar>(
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        Eigen::Vector2d const&,
        DVector2 const&,
        DVector2 const&,
        DVector2 const&,
        DScalar&,
        DScalar&,
        DScalar&);

}
}
