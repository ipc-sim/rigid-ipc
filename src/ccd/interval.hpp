/// An interval object.
#pragma once

#include <boost/numeric/interval.hpp>

namespace ccd {

typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double>>,
        // boost::numeric::interval_lib::checking_no_empty<double>>>
        boost::numeric::interval_lib::checking_base<double>>>
    Interval;

} // namespace ccd
