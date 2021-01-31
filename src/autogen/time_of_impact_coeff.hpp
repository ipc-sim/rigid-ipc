#pragma once

#include <Eigen/Core>
#include <autodiff/autodiff_types.hpp>

namespace ipc::rigid {
namespace autogen {

    template <typename T> using Vector2T = Eigen::Matrix<T, 2, 1>;

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
        T& c);

    template <>
    inline void time_of_impact_coeff<double>(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk,
        const Vector2T<double>& Ui,
        const Vector2T<double>& Uj,
        const Vector2T<double>& Uk,
        double& a,
        double& b,
        double& c)
    {
        a = -Ui[0] * Uj[1] + Ui[0] * Uk[1] + Ui[1] * Uj[0] - Ui[1] * Uk[0]
            - Uj[0] * Uk[1] + Uj[1] * Uk[0];
        b = -Ui[0] * Vj[1] + Ui[0] * Vk[1] + Ui[1] * Vj[0] - Ui[1] * Vk[0]
            + Uj[0] * Vi[1] - Uj[0] * Vk[1] - Uj[1] * Vi[0] + Uj[1] * Vk[0]
            - Uk[0] * Vi[1] + Uk[0] * Vj[1] + Uk[1] * Vi[0] - Uk[1] * Vj[0];
        c = -Vi[0] * Vj[1] + Vi[0] * Vk[1] + Vi[1] * Vj[0] - Vi[1] * Vk[0]
            - Vj[0] * Vk[1] + Vj[1] * Vk[0];

        double max_coeff = std::max(a, std::max(b, c));
        a /= max_coeff;
        b /= max_coeff;
        c /= max_coeff;
    }

} // namespace autogen
} // namespace ipc::rigid

#include "auto_time_of_impact_coeff.ipp"
