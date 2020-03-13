#pragma once

#include <limits>

namespace ccd {
namespace Constants {

    /// \brief Value used to break newton solver ||grad|| < tolerance indicates
    /// we found minimum. Configurable with json file.
    static const double NEWTON_ABSOLUTE_TOLERANCE = 1E-5;

    /// \brief Value used to break linesearch.
    static const double LINESEARCH_MIN_STEP_NORM = 0;

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

    /// \brief ETA used in Etienne Vouga's CCD as an activation distance.
    static const double LINEARIZED_CCD_ETA = 1e-8;

    // ------------------------------------------------------------------------
    // Nonlinear Complementarity Problem
    // ------------------------------------------------------------------------
    /// \brief Number of iterations to take before falling back to normal
    /// direction.
    static const int NCP_FALLBACK_ITERATIONS = 20;
    /// \brief Abolute tolerance on the STIV constraints
    ///        \f$(g(x) ≥ -\epsilon_\text{abs})\f$.
    static const double NCP_ABS_TOL = 1e-12;

    // ------------------------------------------------------------------------
    // Fischer-Newton LCP Solver
    // ------------------------------------------------------------------------
    /// \brief Relative error tolerance for minimizing the Fischer error.
    static const double FISCHER_REL_TOL = 1e-5;
    /// \brief Absolute error tolerance for minimizing the Fischer error.
    static const double FISCHER_ABS_TOL = 1e-28;
    /// \brief Perturbation values used to fix near singular points in
    /// derivative.
    static const double FISCHER_SINGULAR_TOL = 1e-28;
    /// \brief Maximum number of iterations to take while minimizing the Fischer
    /// error.
    static const int FISCHER_MAX_ITER = 3000;

    // ------------------------------------------------------------------------
    // Guass-Seidel LCP Solver
    // ------------------------------------------------------------------------
    /// \brief Maximum number of iterations to take during Guass-Seidel.
    static const int GUASS_SEIDEL_MAX_ITER = 10000;

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
} // namespace ccd
