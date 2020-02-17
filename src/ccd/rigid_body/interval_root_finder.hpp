/// A root finder using interval arithmetic.
#pragma once

#include <functional>

#include <ccd/rigid_body/interval.hpp>

namespace ccd {

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
