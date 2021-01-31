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

        // DO NOT MODIFY THIS FILE. See src/autogen/.tpp for originals.
        const auto helper_0 = Ui[0] * toi - Uj[0] * toi + Vi[0] - Vj[0];
        const auto helper_1 = Ui[1] * toi - Uj[1] * toi + Vi[1] - Vj[1];
        const auto helper_2 = pow(helper_0, 2) + pow(helper_1, 2);
        const auto helper_3 = -helper_0 * (Ui[1] - alpha * (Ui[1] - Uj[1])) + helper_1 * (Ui[0] - alpha * (Ui[0] - Uj[0]));

        assert(sqrt(helper_2) > 0);
        assert(epsilon > 0 || helper_3 != 0);
        volume = (toi - 1.0) * sqrt(pow(epsilon, 2) * helper_2 + pow(helper_3, 2));

        return volume;
    }
}
}