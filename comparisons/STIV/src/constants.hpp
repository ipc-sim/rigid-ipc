#pragma once

#include <limits>

namespace ccd {
namespace Constants {

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

} // namespace Constants
} // namespace ccd
