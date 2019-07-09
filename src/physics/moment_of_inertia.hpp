#pragma once

#include <Eigen/Core>

namespace ccd {
namespace physics {
    ///
    /// \brief moment_of_inertia: computes the 2D moment of intertia wrt to z-axis (screen)
    /// of a set of vertices in body space (i.e centered of mass at 0,0).
    ///
    double moment_of_inertia(
        const Eigen::MatrixXd& vertices, const Eigen::VectorXd& masses);
} // namespace physics
} // namespace ccd
