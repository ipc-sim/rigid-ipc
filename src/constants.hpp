#pragma once

#include <limits>

namespace ipc::rigid {
namespace Constants {

    /// \brief Value used to break newton solver
    static const double DEFAULT_NEWTON_ENERGY_CONVERGENCE_TOL = 1E-5;
    static const double DEFAULT_NEWTON_VELOCITY_CONVERGENCE_TOL = 1E-2;

    /// \brief Value used to break linesearch.
    static const double DEFAULT_LINE_SEARCH_LOWER_BOUND = 1e-12;

    /// \brief Value used to decide when a coeff (a,b,c) is zero.
    static const double TOI_COEFF_EPSILON = 1E-12;

    /// \brief Value used when solving for alpha given a toi, to avoid division
    /// by zero.
    static const double ALPHA_DIVISION_EPSILON = 1E-8;

    /// \brief Distances smaller than this are considered errors.
    static const double POINT_EDGE_SQ_DISTANCE_NEG_TOL = -1E-10;

    /// \brief Value to replace neg-zeros with.
    static const double POINT_EDGE_SQ_DISTANCE_ZERO = 1E-20;

    /// \brief Tolerance when asserting for valid barycentric coordinates.
    static const double PARAMETER_ASSERTION_TOL = 1e-7;

    /// \brief Tolerance of the length of the rigid trajectory in CCD.
    static const double RIGID_CCD_TOI_TOL = 1e-4;
    static const double RIGID_CCD_LENGTH_TOL = 1e-4;

    /// \brief Tolerance on the size of the range of the interval-based CCD.
    static const double INTERVAL_ROOT_FINDER_RANGE_TOL = 1e-8;

    /// \brief Default tolerance used for interval root finding.
    static const double INTERVAL_ROOT_FINDER_TOL = 1e-8;

    /// \brief Default tolerance used for interval root finding.
    static const int INTERVAL_ROOT_FINDER_MAX_ITERATIONS = 10000;

    /// \brief Scaling of κ_min to better condition the system
    static const double DEFAULT_MIN_BARRIER_STIFFNESS_SCALE = 1e11;

    // static const int MAXIMUM_FRICTION_ITERATIONS = 100;

    // ------------------------------------------------------------------------
    // Debugging
    // ------------------------------------------------------------------------
    /// \brief Spacing h of finite differences.
    static const double FINITE_DIFF_H = 1E-7;

    /// \brief Value used for checking if f.d matches.
    static const double FINITE_DIFF_TEST = 1E-4;

    /// \brief Value used for checking if full gradient and assembly of local
    /// gradients match (both computed with autodiff).
    static const double FULL_GRADIENT_TEST = 1E-10;

} // namespace Constants
} // namespace ipc::rigid
