#include "time_of_impact_coeff.hpp"

#include <iostream>

#include <autodiff/autodiff_types.hpp>
#include <autodiff/finitediff.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    T space_time_collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const T& toi, const T& alpha,
        const double epsilon)
    {
        T volume(0.0);

        // {{volume_ccode}}

        return volume;
    }

    template double space_time_collision_volume<double>(
        Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&,
        const double&, const double&, double);

    template DScalar space_time_collision_volume<DScalar>(
        Eigen::Vector2d const&, Eigen::Vector2d const&,
        DVector2 const&, DVector2 const&,
        const DScalar&, const DScalar&, double);

}
}
