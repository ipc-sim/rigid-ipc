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
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const Eigen::VectorX3I& x0,
    Eigen::VectorX3d tol,
    Eigen::VectorX3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& is_domain_valid,
    const Eigen::VectorX3I& x0,
    Eigen::VectorX3d tol,
    Eigen::VectorX3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorX3I(const Eigen::VectorX3I&)>& f,
    const std::function<bool(const Eigen::VectorX3I&)>& constraint_predicate,
    const std::function<bool(const Eigen::VectorX3I&)>& is_domain_valid,
    const Eigen::VectorX3I& x0,
    Eigen::VectorX3d tol,
    Eigen::VectorX3I& x,
    int max_iterations = Constants::INTERVAL_ROOT_FINDER_MAX_ITERATIONS);

} // namespace ipc::rigid
