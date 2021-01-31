#pragma once

#include <autogen/time_of_impact_coeff.hpp>
#include <constants.hpp>

namespace ipc::rigid {

template <typename T>
bool temporal_parameterization_to_spatial(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui,
    const Vector2T<T>& Uj,
    const Vector2T<T>& Uk,
    const T& toi,
    T& alpha)
{
    // solves for alpha
    // Vk + Uk * toi =
    //      (Vi + Ui * toi) + alpha ((Vj + Uj * toi) - (Vi + Ui * toi)).

    auto numerator_0 = Vk[0] - Vi[0] + (Uk[0] - Ui[0]) * toi;
    auto numerator_1 = Vk[1] - Vi[1] + (Uk[1] - Ui[1]) * toi;
    auto denominator_0 = (Vj[0] - Vi[0] + (Uj[0] - Ui[0]) * toi);
    auto denominator_1 = (Vj[1] - Vi[1] + (Uj[1] - Ui[1]) * toi);

    // we need to divide to obtain alpha but we only need one dimension
    // since we want to use autoddif, we use < || > instead of abs
    if (denominator_0 > Constants::ALPHA_DIVISION_EPSILON
        || denominator_0 < -Constants::ALPHA_DIVISION_EPSILON) {
        alpha = numerator_0 / denominator_0;
        return true;
    }
    if (denominator_1 > Constants::ALPHA_DIVISION_EPSILON
        || denominator_1 < -Constants::ALPHA_DIVISION_EPSILON) {
        alpha = numerator_1 / denominator_1;
        return true;
    }
    return false;
}

template <typename T>
bool compute_edge_vertex_time_of_impact(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Vector2T<T>& Ui,
    const Vector2T<T>& Uj,
    const Vector2T<T>& Uk,
    T& toi,
    T& alpha)
{
    T a, b, c;
    autogen::time_of_impact_coeff(Vi, Vj, Vk, Ui, Uj, Uk, a, b, c);

    // At most we will have two solutions for the quadratic equation
    // we initialize them with an invalid value (i.e default failure)
    // our approach uses two different formulas that solve the quadratic
    // equation
    //     (1)         (-b - sqrt(b^2 - 4 * a * c) / (2 * a)
    //     (2)         -2*c / (|b|+ sqrt(b^2 - 4 * a * c))

    auto is_not_zero = [](const T x) -> bool {
        return x > T(Constants::TOI_COEFF_EPSILON)
            || x < -T(Constants::TOI_COEFF_EPSILON);
    };
    T x1(-1), x2(-1);
    T radicand = b * b - 4 * a * c;
    bool a_not_zero = is_not_zero(a);
    bool b_not_zero = is_not_zero(b);
    bool c_not_zero = is_not_zero(c);

    if (radicand > 0) {
        T sqrt_rad = sqrt(radicand);
        if (b > 0) {
            x1 = -2 * c / (b + sqrt_rad);
            // if (a_not_zero) {
            x2 = (-b - sqrt_rad) / (2 * a);
            //}
        } else { // note x1 and x2 switched
            x2 = 2 * c / (-b + sqrt_rad);
            // if (a_not_zero) {
            x1 = (-b + sqrt_rad) / (2 * a);
            //}
        }
    } else if (radicand == 0 && a_not_zero) {
        x1 = (-b) / (2 * a);
    } else if (!a_not_zero && !b_not_zero && !c_not_zero) {
        // check for case a=b=c=0

        // we will use the following approximation:
        //      distance to the closest edge vertex divided by rel
        //      velocity
        Vector2T<T> n = (Vj - Vi).normalized().cast<T>();
        double dist_ik = (Vi - Vk).norm();
        double dist_jk = (Vj - Vk).norm();

        if (dist_ik < dist_jk) { // closest i
            x1 = T(dist_ik) / (Uk - Ui).dot(n);
        } else {
            x1 = T(dist_jk) / (Uk - Uj).dot(-n);
        }
    }
    auto check_solution = [Vi, Vj, Vk, Ui, Uj, Uk](const T& t, T& s) {
        return t >= 0 && t <= 1
            && temporal_parameterization_to_spatial(
                   Vi, Vj, Vk, Ui, Uj, Uk, t, s)
            && s >= 0 && s <= 1;
    };

    // now check which of the solutions are valid
    T alpha1(-1), alpha2(-1);
    bool x1_valid = check_solution(x1, alpha1);
    bool x2_valid = check_solution(x2, alpha2);
    if (x1_valid && x2_valid) {
        toi = std::min(x1, x2);
        alpha = x1 < x2 ? alpha1 : alpha2;
    } else if (x1_valid) {
        toi = x1;
        alpha = alpha1;
    } else {
        toi = x2;
        alpha = alpha2;
    }

    return x1_valid || x2_valid;
}

} // namespace ipc::rigid
