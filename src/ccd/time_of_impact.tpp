#pragma once

#include <autogen/time_of_impact_coeff.hpp>

namespace ccd {
namespace autodiff {

template <typename T>
bool temporal_parameterization_to_spatial(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
    const T& toi, T& alpha)
{
    static const double kEPSILON = 1E-8;

    // solves for alpha
    // Vk + Uk * toi =
    //      (Vi + Ui * toi) + alpha ((Vj + Uj * toi) - (Vi + Ui * toi)).

    auto numerator_0 = Vk[0] - Vi[0] + (Uk[0] - Ui[0]) * toi;
    auto numerator_1 = Vk[1] - Vi[1] + (Uk[1] - Ui[1]) * toi;
    auto denominator_0 = (Vj[0] - Vi[0] + (Uj[0] - Ui[0]) * toi);
    auto denominator_1 = (Vj[1] - Vi[1] + (Uj[1] - Ui[1]) * toi);

    // we need to divide to obtain alpha but we only need one dimension
    // since we want to use autoddif, we use < || > instead of abs
    if (denominator_0 > kEPSILON || denominator_0 < -kEPSILON) {
        alpha = numerator_0 / denominator_0;
        return true;
    }
    if (denominator_1 > kEPSILON || denominator_1 < -kEPSILON) {
        alpha = numerator_1 / denominator_1;
        return true;
    }
    return false;
}

template <typename T>
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
                                        T& toi){
    T alpha;
    return compute_edge_vertex_time_of_impact(Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
}

template <typename T>
bool compute_edge_vertex_time_of_impact(const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
    T& toi, T& alpha)
{
    static const double kEPSILON = 1E-8;

    T a, b, c;
    ccd::autogen::time_of_impact_coeff(Vi, Vj, Vk, Ui, Uj, Uk, a, b, c);

    // At most we will have two solutions for the quadratic equation
    // we initialize them with an invalid value (i.e default failure)
    // our approach uses two different formulas that solve the quadratic
    // equation
    //     (1)         (-b - sqrt(b^2 - 4 * a * c) / (2 * a)
    //     (2)         -2*c / (|b|+ sqrt(b^2 - 4 * a * c))

    T x1(-1), x2(-1);
    T radicand = b * b - 4 * a * c;
    bool a_not_zero = a > kEPSILON || a < -kEPSILON;
    if (radicand > 0) {
        T sqrt_rad = sqrt(radicand);
        if (b > 0) {
            x1 = -2 * c / (b + sqrt_rad);
            if (a_not_zero) {
                x2 = (-b - sqrt_rad) / (2 * a);
            }
        } else { // note x1 and x2 switched
            x2 = 2 * c / (-b + sqrt_rad);
            if (a_not_zero) {
                x1 = (-b + sqrt_rad) / (2 * a);
            }
        }
    }
    // now check which of the solutions are valid
    auto check_solution = [&](T t) {
        return t >= 0 && t <= 1
            && temporal_parameterization_to_spatial(
                Vi, Vj, Vk, Ui, Uj, Uk, t, alpha)
            && alpha >= 0 && alpha <= 1;
    };
    bool x1_valid = check_solution(x1);
    bool x2_valid = check_solution(x2);
    toi = x1_valid && x2_valid ? std::min(x1, x2) : (x1_valid ? x1 : x2);

    return x1_valid || x2_valid;
}

}
}
