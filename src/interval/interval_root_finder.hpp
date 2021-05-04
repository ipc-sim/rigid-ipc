// A root finder using interval arithmetic.
#pragma once

#include <functional>

#include <constants.hpp>
#include <interval/interval.hpp>

namespace ipc::rigid {

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const std::function<bool(const VectorMax3I&)>& is_domain_valid,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<VectorMax3I(const VectorMax3I&)>& f,
    const std::function<bool(const VectorMax3I&)>& constraint_predicate,
    const std::function<bool(const VectorMax3I&)>& is_domain_valid,
    const VectorMax3I& x0,
    VectorMax3d tol,
    VectorMax3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

} // namespace ipc::rigid
