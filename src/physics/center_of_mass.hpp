#pragma once

#include <Eigen/Core>

namespace ccd {
namespace physics {

    Eigen::VectorXd center_of_mass(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges);
}
}
