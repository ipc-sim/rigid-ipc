// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"

#include <stack>

#include <logger.hpp>

namespace ccd {

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x,
    int max_iterations)
{
    return interval_root_finder(
        f, [](const Interval&) { return true; }, x0, tol, x, max_iterations);
}

bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x,
    int max_iterations)
{
    Eigen::VectorX3I x0_vec = Eigen::VectorX3I::Constant(1, x0), x_vec;
    Eigen::VectorX3d tol_vec = Eigen::VectorX3d::Constant(1, tol);
    bool found_root = interval_root_finder(
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return Eigen::VectorX3I::Constant(1, f(x(0)));
        },
        [&](const Eigen::VectorX3I& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        x0_vec, tol_vec, x_vec, max_iterations);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
}

bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x,
    int max_iterations)
{
    return interval_root_finder(
        f, [](const Eigen::VectorX3I&) { return true; }, x0, tol, x,
        max_iterations);
}

inline Eigen::VectorX3d width(const Eigen::VectorX3I& x)
{
    Eigen::VectorX3d w(x.size());
    for (int i = 0; i < x.size(); i++) {
        w(i) = width(x(i));
    }
    return w;
}

inline double diagonal_width(const Eigen::VectorX3I& x)
{
    double w = 0;
    for (int i = 0; i < x.size(); i++) {
        w += width(x(i)) * width(x(i));
    }
    return sqrt(w);
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

void log_octree(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const Eigen::VectorX3I& x0,
    int levels = 5)
{
    if (levels == 1 || !zero_in(f(x0))) {
        std::cout << "[" << logger::fmt_eigen_intervals(x0) << ", "
                  << (zero_in(f(x0)) ? "True" : "False") << "]," << std::endl;
        return;
    }

    for (Interval t : { bisect(x0(0)).first, bisect(x0(0)).second }) {
        for (Interval alpha : { bisect(x0(1)).first, bisect(x0(1)).second }) {
            for (Interval beta :
                 { bisect(x0(2)).first, bisect(x0(2)).second }) {
                Eigen::Vector3I x(t, alpha, beta);
                log_octree(f, x, levels - 1);
            }
        }
    }
}

bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const Eigen::VectorX3I& x0,
    const Eigen::VectorX3d& tol,
    Eigen::VectorX3I& x,
    int max_iterations)
{
    // log_octree(f, x0);

    // Stack of intervals and the last split dimension
    std::stack<std::pair<Eigen::VectorX3I, int>> xs;
    xs.emplace(x0, -1);
    for (int i = 0; i < max_iterations && !xs.empty(); i++) {
        x = xs.top().first;
        int last_split = xs.top().second;
        xs.pop();

        Eigen::VectorX3I y = f(x);

        // spdlog::trace(
        //     "{} ↦ {}", logger::fmt_eigen_intervals(x),
        //     logger::fmt_eigen_intervals(y));

        if (!zero_in(y)) {
            continue;
        }

        Eigen::VectorX3d widths = width(x);
        if ((widths.array() <= tol.array()).all()) {
            if (constraint_predicate(x)) {
                return true;
            }
            continue;
        }

        // Check the diagonal of the range box
        if (diagonal_width(y) <= Constants::INTERVAL_ROOT_FINDER_RANGE_TOL
            && constraint_predicate(x)) {
            return true;
        }

        // Bisect the largest dimension divided by its tolerance
        int split_i = -1;
        for (int i = 0; i < x.size(); i++) {
            if (widths(i) > tol(i)
                && (split_i == -1
                    || widths(i) * tol(split_i) > widths(split_i) * tol(i))) {
                split_i = i;
            }
        }

        std::pair<Interval, Interval> halves = bisect(x(split_i));
        // Push the second half on first so it is examined after the first half
        x(split_i) = halves.second;
        xs.emplace(x, split_i);
        x(split_i) = halves.first;
        xs.emplace(x, split_i);
    }
    if (!xs.empty()) {
        // Return the earliest possible time of impact
        x = xs.top().first;
        xs.pop();
        while (!xs.empty()) {
            if (xs.top().first[0].lower() < x[0].lower()) {
                x = xs.top().first;
            }
            xs.pop();
        }
        return true; // A conservative answer
    }
    return false;
}

} // namespace ccd
