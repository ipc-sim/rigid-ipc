#include "filib_rounding.hpp"

#include <cassert>

#include <filib/fi_lib.h>
#undef local

#include <utils/not_implemented_error.hpp>

namespace ipc::rigid {

double FILibRounding::add_down(double x, double y)
{
    return x == -y ? 0 : q_pred(x + y);
}

double FILibRounding::add_up(double x, double y)
{
    return x == -y ? 0 : q_succ(x + y);
}

double FILibRounding::sub_down(double x, double y)
{
    return x == y ? 0 : q_pred(x - y);
}

double FILibRounding::sub_up(double x, double y)
{
    return x == y ? 0 : q_succ(x - y);
}

double FILibRounding::mul_down(double x, double y)
{
    return (x == 0 || y == 0) ? 0 : q_pred(x * y);
}

double FILibRounding::mul_up(double x, double y)
{
    return (x == 0 || y == 0) ? 0 : q_succ(x * y);
}

double FILibRounding::div_down(double x, double y)
{
    return x == 0 ? 0 : q_pred(x / y);
}

double FILibRounding::div_up(double x, double y)
{
    return x == 0 ? 0 : q_succ(x / y);
}

double FILibRounding::sqrt_down(double x)
{
    return x == 0 ? 0 : r_pred(q_sqrt(x));
}

double FILibRounding::sqrt_up(double x)
{
    return x == 0 ? 0 : r_succ(q_sqrt(x));
}

double FILibRounding::exp_down(double x)
{
    double r = x == 0 ? 1.0 : (x <= q_mine ? 0 : (q_exp(x) * q_exem));
    // Snap negative r values to 0
    if (r < 0) {
        r = 0;
    }
    // Snap negative r values to 1 if x is positive
    if (x >= 0 && r < 1.0) {
        r = 1;
    }
    return r;
}

double FILibRounding::exp_up(double x)
{
    double r = x == 0 ? 1.0 : (x <= q_mine ? q_minr : (q_exp(x) * q_exep));
    // Snap negative r values to 0
    if (r < 0) {
        r = 0;
    }
    // Snap negative r values to 1 if x is positive
    if (x >= 0 && r < 1.0) {
        r = 1;
    }
    return r;
}

// double FILibRounding::log_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("");
// }

// double FILibRounding::log_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("");
// }

double FILibRounding::cos_down(double x)
{
    double r;
    if (x < -q_sint[2] || x > q_sint[2]) {
        r = -1.0;
    } else {
        r = q_cos(x);
        r *= r < 0 ? q_cosp : q_cosm;
    }
    if (r < -1.0) {
        r = -1;
    }
    if (r > 1.0) {
        r = 1;
    }
    return r;
}

double FILibRounding::cos_up(double x)
{
    double r;
    if (x < -q_sint[2] || x > q_sint[2]) {
        r = 1.0;
    } else {
        r = q_cos(x);
        r *= r < 0 ? q_cosm : q_cosp;
    }
    if (r < -1.0) {
        r = -1;
    }
    if (r > 1.0) {
        r = 1;
    }
    return r;
}

// double FILibRounding::tan_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("tan_down is not implemented!");
// }

// double FILibRounding::tan_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("tan_up is not implemented!");
// }

// double FILibRounding::asin_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("asin_down is not implemented!");
// }

// double FILibRounding::asin_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("asin_up is not implemented!");
// }

// double FILibRounding::acos_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("acos_down is not implemented!");
// }

// double FILibRounding::acos_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("acos_up is not implemented!");
// }

// double FILibRounding::atan_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("atan_down is not implemented!");
// }

// double FILibRounding::atan_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("atan_up is not implemented!");
// }

// double FILibRounding::sinh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("sinh_down is not implemented!");
// }

// double FILibRounding::sinh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("sinh_up is not implemented!");
// }

// double FILibRounding::cosh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("cosh_down is not implemented!");
// }

// double FILibRounding::cosh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("cosh_up is not implemented!");
// }

// double FILibRounding::tanh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("tanh_down is not implemented!");
// }

// double FILibRounding::tanh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("tanh_up is not implemented!");
// }

// double FILibRounding::asinh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("asinh_down is not implemented!");
// }

// double FILibRounding::asinh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("asinh_up is not implemented!");
// }

// double FILibRounding::acosh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("acosh_down is not implemented!");
// }

// double FILibRounding::acosh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("acosh_up is not implemented!");
// }

// double FILibRounding::atanh_down(double x)
// {
//     assert(false);
//     throw NotImplementedError("atanh_down is not implemented!");
// }

// double FILibRounding::atanh_up(double x)
// {
//     assert(false);
//     throw NotImplementedError("atanh_up is not implemented!");
// }

double FILibRounding::median(double x, double y)
{
    this->to_nearest();
    return this->force_rounding(q_mid({ x, y }));
}

double FILibRounding::int_down(double x)
{
    this->downward();
    return this->to_int(x);
}

double FILibRounding::int_up(double x)
{
    this->upward();
    return this->to_int(x);
}

} // namespace ipc::rigid
