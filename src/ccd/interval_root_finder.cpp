// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"

#include <logger.hpp>

namespace ccd {

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x)
{
    return interval_root_finder(
        f, [](const Interval&) { return true; }, x0, tol, x);
}

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x)
{
    Eigen::VectorXI x0_vec = Eigen::VectorXI::Constant(1, x0), x_vec;
    Eigen::VectorXd tol_vec = Eigen::VectorXd::Constant(1, tol);
    bool found_root = interval_root_finder(
        [&](const Eigen::VectorXI& x) {
            assert(x.size() == 1);
            return Eigen::VectorXI::Constant(1, f(x(0)));
        },
        [&](const Eigen::VectorXI& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        x0_vec, tol_vec, x_vec);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
}

bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const Eigen::VectorXI& x0,
    const Eigen::VectorXd& tol,
    Eigen::VectorXI& x)
{
    return interval_root_finder(
        f, [](const Eigen::VectorXI&) { return true; }, x0, tol, x);
}

Eigen::ArrayXd width(const Eigen::VectorXI& x)
{
    Eigen::ArrayXd w(x.size());
    for (int i = 0; i < x.size(); i++) {
        w(i) = width(x(i));
    }
    return w;
}

template <int dim, int max_dim = dim>
inline bool zero_in(Eigen::Vector<Interval, dim, max_dim> X)
{
    // Check if the origin is in the n-dimensional interval
    for (int i = 0; i < X.size(); i++) {
        if (!boost::numeric::zero_in(X(i))) {
            return false;
        }
    }
    return true;
}

bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const std::function<bool(const Eigen::VectorXI&)>& constraint_predicate,
    const Eigen::VectorXI& x0,
    const Eigen::VectorXd& tol,
    Eigen::VectorXI& x)
{
    // Stack of intervals and the last split dimension
    std::stack<std::pair<Eigen::VectorXI, int>> xs;
    xs.emplace(x0, -1);
    while (!xs.empty()) {
        x = xs.top().first;
        int last_split = xs.top().second;
        xs.pop();

        Eigen::VectorXI y = f(x);

        if (!zero_in(y)) {
            continue;
        }

        Eigen::ArrayXd widths = width(x);
        if ((widths <= tol.array()).all()) {
            if (constraint_predicate(x)) {
                return true;
            }
            continue;
        }

        // Bisect the next dimension that is greater than its tolerance
        int split_i;
        for (int i = 1; i <= x.size(); i++) {
            split_i = (last_split + i) % x.size();
            if (widths(split_i) > tol(split_i)) {
                break;
            }
        }
        std::pair<Interval, Interval> halves = bisect(x(split_i));
        Eigen::VectorXI x1 = x;
        // Push the second half on first so it is examined after the first half
        x(split_i) = halves.second;
        xs.emplace(x, split_i);
        x(split_i) = halves.first;
        xs.emplace(x, split_i);
    }
    return false;
}

} // namespace ccd
