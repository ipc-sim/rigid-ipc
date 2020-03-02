/// An interval object.
#pragma once

#include <boost/numeric/interval.hpp>

namespace ccd {

namespace interval_options {
    typedef boost::numeric::interval_lib::checking_base<double> CheckingPolicy;
} // namespace interval_options

#if defined(__APPLE__)

// clang-format off
#warning "Rounding modes seem to be broken with trigonometric functions on macOS, unable to compute exact interval arithmetic!"
// clang-format on
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_exact<double>>,
        interval_options::CheckingPolicy>>
    Interval;

#else

// Use proper rounding arithmetic
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double>>,
        interval_options::CheckingPolicy>>
    Interval;

#endif

} // namespace ccd
