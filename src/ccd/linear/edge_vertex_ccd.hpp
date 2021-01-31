#pragma once

#include <Eigen/Core>
#include <autodiff/autodiff_types.hpp>

/**
 * @namespace ipc::rigid:
 * @brief
 */
namespace ipc::rigid {

template <typename T> using Vector2T = Eigen::Matrix<T, 2, 1>;

/**
 * @brief Given toi returns the position \f$\alpha\f$ along the
 * \f$edge_{ij}\f$ where the impact between vertex \f$k\f$ and
 * \f$edge_{ij}\f$  takes place.
 *
 *   @param[in]   V_{i,j,k}       : vertex positions
 *   @param[in]   U_{i,j,k}       : vertex displacements
 *   @param[in]   toi             : time of impact
 *
 *   @param[out]  alpha          : position along the \f$edge_{ij}\f$ where
 * the impact takes place
 *
 * @return true if a valid \f$\alpha\f$ exists.
 */
template <typename T>
bool temporal_parameterization_to_spatial(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui,
    const Vector2T<T>& Uj,
    const Vector2T<T>& Uk,
    const T& toi,
    T& alpha);

/**
 * Computes the time of impact between an \f$edge_{ij}\f$ and vertex \f$k\f$
 *
 *   @param[in]   V_{i,j,k}       : vertex positions
 *   @param[in]   U_{i,j,k}       : vertex displacements
 *
 *   @param[out]  toi             : time of FIRST impact
 *
 * @return true if an impact happens at time \f$ t \in [0, 1]\f$
 */
template <typename T>
bool compute_edge_vertex_time_of_impact(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui,
    const Vector2T<T>& Uj,
    const Vector2T<T>& Uk,
    T& toi,
    T& alpha);

template <typename T>
bool compute_edge_vertex_time_of_impact(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui,
    const Vector2T<T>& Uj,
    const Vector2T<T>& Uk,
    T& toi)
{
    T alpha;
    return compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
}

} // namespace ipc::rigid

#include "edge_vertex_ccd.tpp"
