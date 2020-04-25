// A root finder using interval arithmetic.
#pragma once

#include <functional>

#include <ccd/interval.hpp>
#include <constants.hpp>

namespace ccd {

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const Interval& x0,
    double tol,
    Interval& x);

/// Find the first root of a function f: I ↦ I
bool interval_root_finder(
    const std::function<Interval(const Interval&)>& f,
    const std::function<bool(const Interval&)>& constraint_predicate,
    const Interval& x0,
    double tol,
    Interval& x);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const Eigen::VectorXI& x0,
    const Eigen::VectorXd& tol,
    Eigen::VectorXI& x);

/// Find if the origin is in the range of a function f: Iⁿ ↦ Iⁿ
bool interval_root_finder(
    const std::function<Eigen::VectorXI(const Eigen::VectorXI&)>& f,
    const std::function<bool(const Eigen::VectorXI&)>& constraint_predicate,
    const Eigen::VectorXI& x0,
    const Eigen::VectorXd& tol,
    Eigen::VectorXI& x);

} // namespace ccd
