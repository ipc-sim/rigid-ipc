#pragma once

#include <boost/numeric/interval.hpp>
#include <functional>

namespace ccd {

typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double>>,
        boost::numeric::interval_lib::checking_base<double>>>
    Interval;

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    Interval& x,
    double tol = 1e-12);

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    Interval& x,
    double tol = 1e-12);

} // namespace ccd
