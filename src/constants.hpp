#pragma once

namespace ccd {
namespace Constants {


    /// \brief NEWTON_ABSOLUTE_TOLERANCE used to break newton solver
    /// ||grad|| < tolerance indicates we found minimum. Configurable
    /// with json file.
    static const double NEWTON_ABSOLUTE_TOLERANCE = 1E-5;

    /// \brief NEWTON_MIN_STEP_NORM used to break linesearch
    static const double LINESEARCH_MIN_STEP_NORM = 0;

    // Used for debugging
    // --------------------------------------------------------
    /// \brief FINITE_DIFF_H spacing h of finite differences
    static const double FINITE_DIFF_H = 1E-7;

    /// \brief FINITE_DIFF_TEST used for checking if f.d matches
    static const double FINITE_DIFF_TEST = 1E-4;

    /// \brief FULL_GRADIENT_TEST used for checking if full gradient
    /// and assembly of local gradients match (both computed with autodiff)
    static const double FULL_GRADIENT_TEST = 1E-16;

}

} // namespace ccd
