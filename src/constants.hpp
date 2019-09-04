#pragma once

namespace ccd {
namespace Constants {


    /// \brief NEWTON_ABSOLUTE_TOLERANCE used to break newton solver
    /// ||grad|| < tolerance indicates we found minimum. Configurable
    /// with json file.
    static const double NEWTON_ABSOLUTE_TOLERANCE = 1E-5;

    /// \brief NEWTON_MIN_STEP_NORM used to break linesearch
    static const double LINESEARCH_MIN_STEP_NORM = 0;

    /// \brief TOI_COEFF_EPSILON used to decide when a coeff (a,b,c) is zero.
    static const double TOI_COEFF_EPSILON = 1E-8;

    /// \brief ALPHA_DIVISION_EPSILON when solving for alpha given a toi,
    /// to avoid division by zero
    static const double ALPHA_DIVISION_EPSILON = 1E-8;

    /// \brief POINT_EDGE_SQ_DISTANCE_NEG_TOL distances smaller than this
    /// are considered errors
    static const double POINT_EDGE_SQ_DISTANCE_NEG_TOL = -1E-10;

    /// \brief POINT_EDGE_SQ_DISTANCE_ZERO value to replace neg-zeros
    /// with
    static const double POINT_EDGE_SQ_DISTANCE_ZERO = 1E-20;

    // Used for debugging
    // --------------------------------------------------------
    /// \brief FINITE_DIFF_H spacing h of finite differences
    static const double FINITE_DIFF_H = 1E-7;

    /// \brief FINITE_DIFF_TEST used for checking if f.d matches
    static const double FINITE_DIFF_TEST = 1E-4;

    /// \brief FULL_GRADIENT_TEST used for checking if full gradient
    /// and assembly of local gradients match (both computed with autodiff)
    static const double FULL_GRADIENT_TEST = 1E-10;

}

} // namespace ccd
