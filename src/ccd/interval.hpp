/// An interval object.
#pragma once

#include <boost/numeric/interval.hpp>

namespace ccd {

namespace interval_options {
    // typedef boost::numeric::interval_lib::checking_no_empty<double>
    //     CheckingPolicy;
    typedef boost::numeric::interval_lib::checking_base<double> CheckingPolicy;
} // namespace interval_options

#if defined(__clang__)

// clang-format off
#warning "Exact interval arithmetic with rounding mode is not supported with clang!"
// clang-format on
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_exact<double>>,
        interval_options::CheckingPolicy>>
    Interval;

#elif defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1310)

// Use proper rounding arithmetic
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double>>,
        interval_options::CheckingPolicy>>
    Interval;

#else

typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_exact<double>>,
        interval_options::CheckingPolicy>>
    Interval;

#endif

} // namespace ccd
