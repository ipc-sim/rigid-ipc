// Test Pose

#include <catch2/catch.hpp>

#include <physics/pose.hpp>

TEST_CASE("Poses to dofs", "[physics][pose]")
{
    using namespace ccd::physics;
    int dim = GENERATE(2, 3);
    int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);
    Eigen::VectorXd dofs =
        Eigen::VectorXd(num_bodies * Pose<double>::dim_to_ndof(dim));
    std::vector<Pose<double>> poses = Pose<double>::dofs_to_poses(dofs, dim);
    Eigen::VectorXd returned_dofs = Pose<double>::poses_to_dofs(poses);
    CHECK((dofs - returned_dofs).squaredNorm() == Approx(0));
}
