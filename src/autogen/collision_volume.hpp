#pragma once

#include <Eigen/Core>
#include <autodiff/autodiff_types.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    using Vector2T = Eigen::Matrix<T, 2, 1>;

    /**
    * Compute the volume of intersection for an edge given a time
    * of intersection (toi) and position of intersection (alpha)
    *
    * \f$V = (1-\tau_I)\sqrt{\epsilon^2 \|e(\tau_I)\|^2 + (U_{ij} \cdot
    * e(\tau_I)^\perp)^2}\f$
    *
    *     @param V_{ijk}          : Vertices positions.
    *     @param U_{ijk}          : Vertices displacements.
    *     @param epsilon          : The time scale used for minimal volume.
    *
    * @return                     : The space-time interference volume.
    */
    template <typename T>
    T space_time_collision_volume(
        const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Vector2T<T>& Ui,
        const Vector2T<T>& Uj,
        const T& toi, const T& alpha,
        const double epsilon);

}
}

#include "auto_collision_volume.ipp"
