#pragma once

#include "time_of_impact_coeff.hpp"


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
}
}
