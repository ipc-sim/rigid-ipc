// A root finder using interval arithmetic.
#include "interval_root_finder.hpp"

#include <stack>

#include <logger.hpp>

namespace ipc::rigid {

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
    VectorMax3I x0_vec = VectorMax3I::Constant(1, x0), x_vec;
    VectorMax3d tol_vec = VectorMax3d::Constant(1, tol);
    bool found_root = interval_root_finder(
        [&](const VectorMax3I& x) {
            assert(x.size() == 1);
            return VectorMax3I::Constant(1, f(x(0)));
        },
        [&](const VectorMax3I& x) {
            assert(x.size() == 1);
            return constraint_predicate(x(0));
        },
        [&](const VectorMax3I& x) { return true; }, x0_vec, tol_vec, x_vec,
        max_iterations);
    if (found_root) {
        assert(x_vec.size() == 1);
        x = x_vec(0);
    }
    return found_root;
}

bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations)
{
    return interval_root_finder(
        f, [](const VectorMax3I&) { return true; }, x0, tol, x, max_iterations);
}

bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const std::function<bool(const VectorMax3I&)>& is_domain_valid,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations)
{
    return interval_root_finder(
        f, [](const VectorMax3I&) { return true; }, is_domain_valid, x0, tol, x,
        max_iterations);
}

void log_octree(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const VectorMax3I& x0,
    int levels = 5)
{
    if (levels == 1 || !zero_in(f(x0))) {
        std::cout << "[" << fmt_eigen_intervals(x0) << ", "
                  << (zero_in(f(x0)) ? "True" : "False") << "]," << std::endl;
        return;
    }

    std::pair<Interval, Interval> t_split = bisect(x0(0));
    for (Interval t : { t_split.first, t_split.second }) {
        std::pair<Interval, Interval> alpha_split = bisect(x0(1));
        for (Interval alpha : { alpha_split.first, alpha_split.second }) {
            std::pair<Interval, Interval> beta_split = bisect(x0(2));
            for (Interval beta : { beta_split.first, beta_split.second }) {
                Vector3I x(t, alpha, beta);
                log_octree(f, x, levels - 1);
            }
        }
    }
}

bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const std::function<bool(const VectorMax3I&)>& constraint_predicate,
    const std::function<bool(const VectorMax3I&)>& is_domain_valid,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations)
{
    // log_octree(f, x0);

    // Keep searching for earlier roots (assumes time is first coordinate)
    VectorMax3I earliest_root = VectorMax3I::Constant(
        x0.size(), Interval(std::numeric_limits<double>::infinity()));
    bool found_root = false;

    // Stack of intervals and the last split dimension
    std::stack<VectorMax3I> xs;
    xs.push(x0);

    // If the start is a root then we are in trouble, so we should reduce the
    // tolerance.
    VectorMax3I x_tol(tol.size());
    for (int i = 0; i < x_tol.size(); i++) {
        x_tol(i) = Interval(0, tol(i));
    }
    if (zero_in(f(x_tol))) {
        tol(0) /= 1e2;
    }

    // TODO: Enable max_iterations
    for (size_t iter = 0; !xs.empty(); iter++) {
        x = xs.top();
        xs.pop();

        // Skip any interval that is not before the earliest root
        if (x[0].lower() >= earliest_root[0].lower()) {
            continue;
        }

        if (!is_domain_valid(x)) {
            continue;
        }

        VectorMax3I y = f(x);

        // spdlog::critical(
        //     "{} ↦ {}", fmt_eigen_intervals(x),
        //     fmt_eigen_intervals(y));

        if (!zero_in(y)) {
            continue;
        }

        VectorMax3d widths = width(x);
        bool all_tol_sat = (widths.array() <= tol.array()).all();
        bool all_widths_zero = (widths.array() <= 1e-10).all();
        if ((x[0].lower() > 0 || all_widths_zero) && all_tol_sat) {
            if (constraint_predicate(x)) {
                earliest_root = x;
                found_root = true;
            }
            continue;
        }

        // Check the diagonal of the range box
        // if (diagonal_width(y) <= Constants::INTERVAL_ROOT_FINDER_RANGE_TOL
        //     && constraint_predicate(x)) {
        //     earliest_root = x;
        //     found_root = true;
        //     continue;
        // }

        // Bisect the largest dimension divided by its tolerance
        int split_i = -1;
        for (int i = 0; i < x.size(); i++) {
            if ((all_tol_sat || widths(i) > tol(i))
                && (split_i == -1
                    || widths(i) * tol(split_i) > widths(split_i) * tol(i))) {
                split_i = i;
            }
        }
        assert(split_i >= 0 && split_i <= x.size());

        std::pair<Interval, Interval> halves = bisect(x(split_i));
        // Push the second half on first so it is examined after the first half
        x(split_i) = halves.second;
        xs.push(x);
        x(split_i) = halves.first;
        xs.push(x);
    }
    // TODO: This is needs to be updated when max_iterations is renabled
    // if (!xs.empty()) {
    //     // Return the earliest possible time of impact
    //     x = xs.top().first;
    //     xs.pop();
    //     while (!xs.empty()) {
    //         if (xs.top().first[0].lower() < x[0].lower()) {
    //             x = xs.top().first;
    //         }
    //         xs.pop();
    //     }
    //     return true; // A conservative answer
    // }
    x = earliest_root;
    return found_root;
}

} // namespace ipc::rigid
