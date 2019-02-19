#include "collision_volume.hpp"

#include <iostream>

#include <autodiff/finitediff.hpp>

namespace ccd {
namespace autogen {

    template <typename T>
    bool temporal_parameterization_to_spatial(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const T& toi, T& alpha)
    {
        static const double kEPSILON = 1E-8;

        // solves for alpha
        // Vk + Uk * toi = (Vi + Ui * toi) + alpha ((Vj + Uj * toi) - (Vi + Ui *
        // toi)).
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
        T& toi)
    {
        static const double kEPSILON = 1E-8;
        // We want to solve the equation
        //      a * (toi)^2 + b * (toi) + c = 0
        // the expressions for the coefficients are filled in using sympy

        T a, b, c;
        // {{toi_abc_ccode}}

        T alpha;
        auto check_solution = [&](T t) {
            return t >= 0 && t <= 1
                && temporal_parameterization_to_spatial(
                       Vi, Vj, Vk, Ui, Uj, Uk, t, alpha)
                && alpha >= 0 && alpha <= 1;
        };

        // Solve the quadratic equation
        if (a > kEPSILON
            || a < -kEPSILON) { //         the equation is quadratic
            // a*t^2 + b*t + c = 0 => t = (-b +/- sqrt(b^2 - 4ac)) / 2a
            T radicand = b * b - 4 * a * c;
            if (radicand < 0) {
                return false;
            }
            T sqrt_rad = sqrt(radicand);

            // the two candidates
            T toi_neg = (-b - sqrt_rad) / (2 * a);
            T toi_pos = (-b + sqrt_rad) / (2 * a);

            bool toi_neg_valid = check_solution(toi_neg);
            bool toi_pos_valid = check_solution(toi_pos);
            if (toi_neg_valid && toi_pos_valid) {
                toi = std::min(toi_neg, toi_pos);
                return true;
            } else if (toi_neg_valid) {
                toi = toi_neg;
                return true;
            } else {
                toi = toi_pos;
                return toi_pos_valid;
            }

        } else if (b > kEPSILON || b < -kEPSILON) { //  the equation is linear
            // b*t + c = 0 => t = -c / b
            toi = -c / b;
            return check_solution(toi);

        } else if (c > kEPSILON
            || c < -kEPSILON) { //  No solution (c !=0 && c == 0)
            return false;

        } else { //                                     Infinite solutions

            // This case happens if the trajectories of Vi,j,k are in the same
            // line as the edge For an initial, contact free state, the position
            // of Vk along the edge must be either alpha < 0 or > alpha 1
            T t0(0.0), alpha_t0;
            temporal_parameterization_to_spatial(
                Vi, Vj, Vk, Ui, Uj, Uk, t0, alpha_t0);
            // For a contact to ocurr, the final position might be `oppossite`
            // the initial position
            T t1(1.0), alpha_t1;
            temporal_parameterization_to_spatial(
                Vi, Vj, Vk, Ui, Uj, Uk, t1, alpha_t1);

            auto alpha0 = alpha_t0 < 0.0
                && alpha_t1 >= 0.0; // first contact is at alpha = 0.0
            auto alpha1 = alpha_t0 > 1.0
                && alpha_t1 <= 1.0; // first contact is at alpha = 1.0
            if (alpha0 || alpha1) {

                alpha = alpha0 ? 0.0 : 1.0;
                T delta_alpha = alpha_t1 - alpha_t0;
                auto delta_non_zero
                    = delta_alpha > kEPSILON || delta_alpha < kEPSILON;
                toi = delta_non_zero ? (alpha - alpha_t0) / delta_alpha
                                     : T(0.0);
                return true;

            } else if (alpha_t0 > 0.0 && alpha_t0 < 1.0) {
                // initial state is not contact free
                // we set first contact is at time 0
                toi = 0.0;
                alpha = alpha_t0;
                return true;
            } else {
                return false;
            }
        }
    }

    template <typename T>
    bool compute_edge_vertex_time_of_impact_v2(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        T& toi)
    {
        static const double kEPSILON = 1E-8;
        // We want to solve the equation
        //      a * (toi)^2 + b * (toi) + c = 0
        // the expressions for the coefficients are filled in using sympy

        T a, b, c;
        // {{toi_abc_ccode}}

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

        T alpha;
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

    bool compute_edge_vertex_time_of_impact_grad(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk, const int version, Vector8d& grad)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DScalar(6, 0.0);
        DScalar(7, 0.0);

        DScalar toi(0.0);

        bool success = false;
        if (version == 1) {
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vk, DUi, DUj, DUk, toi);
        } else {
            success = compute_edge_vertex_time_of_impact_v2(
                Vi, Vj, Vk, DUi, DUj, DUk, toi);
        }
        if (success) {
            grad = toi.getGradient();
        }
        return success;
    }

    bool compute_edge_vertex_time_of_impact_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Ui, const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& Uk, Vector8d& grad)
    {
        auto f = [Vi, Vj, Vk](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);

            double toi;
            bool success = ccd::autogen::compute_edge_vertex_time_of_impact_v2(
                Vi, Vj, Vk, ui, uj, uk, toi);
            return success ? toi : 0.0;
        };

        Eigen::VectorXd finite_diff;
        Vector8d x;
        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) << 0.0, 0.0;
        ccd::finite_gradient(x, f, finite_diff);
        grad << finite_diff;
        return true;
    }

    template <typename T>
    T _collision_volume_(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const double epsilon)
    {
        using namespace std;
        double alpha = 0.5, toi = 0.5;
        T volume(0.0);

        // {{volume_ccode}}

        return volume;
    }

    double collision_volume(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const double epsilon)
    {
        return _collision_volume_(Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, epsilon);
    }

    Vector8d collision_volume_grad(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const double epsilon)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DVector2 DUl = dvector(6, Ul);

        DScalar volume(0.0);

        volume
            = _collision_volume_(Vi, Vj, Vk, Vl, DUi, DUj, DUk, DUl, epsilon);

        return volume.getGradient();
    }

    template bool compute_edge_vertex_time_of_impact<double>(
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, double&);

    template bool compute_edge_vertex_time_of_impact_v2<double>(
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&,
        Eigen::Matrix<double, 2, 1, 0, 2, 1> const&, double&);
}
}
