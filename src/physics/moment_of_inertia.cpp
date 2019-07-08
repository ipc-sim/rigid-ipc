#include "moment_of_inertia.hpp"

namespace ccd {
namespace physics {
    double moment_of_inertia(
        const Eigen::MatrixXd& vertices, const Eigen::VectorXd& masses)
    {
        // NOTE: this assumes center of mass is at 0,0
        // NOTE: this only works in 2D as we take the moment of intertia wrt z-axis
        //
        // \Sum m_i * r_i \dot r_i
        return (masses.array() * vertices.rowwise().squaredNorm().array()).sum();
    }
} // namespace physics
} // namespace ccd
