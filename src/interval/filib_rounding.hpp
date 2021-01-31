#pragma once

#include <boost/numeric/interval.hpp>

namespace ipc::rigid {

// A wrapper for filibs rounding
struct FILibRounding : boost::numeric::interval_lib::rounding_control<double> {
    // default constructor, destructor
    // FILibRounding() {}
    // ~FILibRounding() {}
    void init() {}
    // mathematical operations
    double add_down(double x, double y); // [-∞;+∞][-∞;+∞]
    double add_up(double x, double y);   // [-∞;+∞][-∞;+∞]
    double sub_down(double x, double y); // [-∞;+∞][-∞;+∞]
    double sub_up(double x, double y);   // [-∞;+∞][-∞;+∞]
    double mul_down(double x, double y); // [-∞;+∞][-∞;+∞]
    double mul_up(double x, double y);   // [-∞;+∞][-∞;+∞]
    double div_down(double x, double y); // [-∞;+∞]([-∞;+∞]-{0})
    double div_up(double x, double y);   // [-∞;+∞]([-∞;+∞]-{0})
    double sqrt_down(double x);          // ]0;+∞]
    double sqrt_up(double x);            // ]0;+∞]
    double exp_down(double x);           // [-∞;+∞]
    double exp_up(double x);             // [-∞;+∞]
    // double log_down(double x);           // ]0;+∞]
    // double log_up(double x);             // ]0;+∞]
    double cos_down(double x); // [0;2π]
    double cos_up(double x);   // [0;2π]
    // double tan_down(double x);           // ]-π/2;π/2[
    // double tan_up(double x);             // ]-π/2;π/2[
    // double asin_down(double x);          // [-1;1]
    // double asin_up(double x);            // [-1;1]
    // double acos_down(double x);          // [-1;1]
    // double acos_up(double x);            // [-1;1]
    // double atan_down(double x);          // [-∞;+∞]
    // double atan_up(double x);            // [-∞;+∞]
    // double sinh_down(double x);          // [-∞;+∞]
    // double sinh_up(double x);            // [-∞;+∞]
    // double cosh_down(double x);          // [-∞;+∞]
    // double cosh_up(double x);            // [-∞;+∞]
    // double tanh_down(double x);          // [-∞;+∞]
    // double tanh_up(double x);            // [-∞;+∞]
    // double asinh_down(double x);         // [-∞;+∞]
    // double asinh_up(double x);           // [-∞;+∞]
    // double acosh_down(double x);         // [1;+∞]
    // double acosh_up(double x);           // [1;+∞]
    // double atanh_down(double x);         // [-1;1]
    // double atanh_up(double x);           // [-1;1]

    /// @brief median of two numbers rounded to closest number
    double median(double x, double y); // [-∞;+∞][-∞;+∞]

    /// @brief round down to integer
    double int_down(double x); // [-∞;+∞]
    /// @brief round up to integer
    double int_up(double x); // [-∞;+∞]

    /// @brief convert to the closest double (rounding down)
    template <typename T> double conv_down(const T& x)
    {
        this->downward();
        return this->force_rounding(x);
    }
    /// @brief convert to the closest double (rounding up)
    template <typename T> double conv_up(const T& x)
    {
        this->upward();
        return this->force_rounding(x);
    }
};

} // namespace ipc::rigid
