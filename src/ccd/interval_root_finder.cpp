// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"

#include <logger.hpp>

namespace ccd {

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    Interval& x,
    double tol)
{
    return interval_root_finder(
        f, [](const Interval&) { return true; }, x0, x, tol);
}

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    Interval& x,
    double tol)
{
    Eigen::VectorXI x0_vec = Eigen::VectorXI::Constant(1, x0), x_vec;
    bool found_root = interval_root_finder(
        [&](const Eigen::VectorXI& x) {
            assert(x.size() == 1);
            return Eigen::VectorXI::Constant(1, f(x(0)));
        },
        [&](const Eigen::VectorXI& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        x0_vec, x_vec, tol);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
}

double width(const Eigen::VectorXI& x)
{
    double w = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < x.size(); i++) {
        w = std::max(w, width(x(i)));
    }
    return w;
}

bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const Eigen::VectorXI& x0,
    Eigen::VectorXI& x,
    double tol,
    int last_split)
{
    return interval_root_finder(
        f, [](const Eigen::VectorXI&) { return true; }, x0, x, tol, last_split);
}

bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const std::function<bool(const Eigen::VectorXI&)>& constraint_predicate,
    const Eigen::VectorXI& x0,
    Eigen::VectorXI& x,
    double tol,
    int last_split)
{
    Eigen::VectorXI y = f(x0);
    assert(y.size() == x0.size());
    // spdlog::debug(
    //     "{} ↦ {}", logger::fmt_eigen_intervals(x0),
    //     logger::fmt_eigen_intervals(y));
    for (int i = 0; i < y.size(); i++) {
        // Check that every interval contains zero
        if (!zero_in(y(i))) {
            return false;
        }
    }
    if (width(x0) <= tol) {
        x = x0;
        return constraint_predicate(x);
    }

    // Bisect the dimension
    int split_i = (last_split + 1) % x0.size();
    std::pair<Interval, Interval> halves = bisect(x0(split_i));
    Eigen::VectorXI x1 = x0;
    x1(split_i) = halves.first;
    // Check the first half of the box
    if (!interval_root_finder(f, constraint_predicate, x1, x, tol, split_i)) {
        x1(split_i) = halves.second;
        // Check the second half of the box
        return interval_root_finder(
            f, constraint_predicate, x1, x, tol, split_i);
    }
    return true;
}

} // namespace ccd
