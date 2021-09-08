#pragma once

#include "time_of_impact_coeff.hpp"

#include <iostream>

#include <autodiff/autodiff_types.hpp>

namespace ipc::rigid {
namespace autogen {

    template <typename T>
    void time_of_impact_coeff(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const Vector2T<T>& Uk,
        T& a,
        T& b,
        T& c)
    {
        // {{toi_abc_ccode}}
    }

} // namespace autogen
} // namespace ipc::rigid
