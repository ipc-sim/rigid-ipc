#include "interval_root_finder.hpp"

#include <iostream>

namespace ccd {

bool interval_root_finder(const std::function<Interval(Interval)>& f,
    const Interval& x0,
    Interval& x,
    double tol)
{
    // std::cout << "[" << x0.lower() << ", " << x0.upper() << "] ↦ ";
    Interval y = f(x0);
    // std::cout << "[" << y.lower() << ", " << y.upper() << "]" << std::endl;
    if (!zero_in(y)) {
        x = Interval::empty(); // Return an empty interval
        return false;
    }
    if (width(x0) <= tol) {
        x = x0;
        return true;
    }
    std::pair<Interval, Interval> halves = bisect(x0);
    return interval_root_finder(f, halves.first, x, tol)
        || interval_root_finder(f, halves.second, x, tol);
}

} // namespace ccd
