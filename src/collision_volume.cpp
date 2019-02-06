#include <FixingCollisions/collision_volume.hpp>


namespace ccd {

double collision_volume(
    const Eigen::MatrixX2d& /*vertices*/,
    const Eigen::MatrixX2d& /*displacements*/,
    const Eigen::MatrixX2i& /*edges*/,
    const Impact& /*impact*/,
    const Eigen::VectorXd /*grad_volume*/)
{
    return 0.0;
}
}
