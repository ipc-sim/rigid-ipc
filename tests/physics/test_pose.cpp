// Test Pose

#include <catch2/catch.hpp>

#include <physics/pose.hpp>

TEST_CASE("Poses to dofs", "[physics][pose]")
{
    using namespace ccd::physics;
    int dim = GENERATE(2, 3);
    int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);
    Eigen::VectorXd dofs =
        Eigen::VectorXd::Random(num_bodies * Pose<double>::dim_to_ndof(dim));
    Poses<double> poses = Pose<double>::dofs_to_poses(dofs, dim);
    Eigen::VectorXd returned_dofs = Pose<double>::poses_to_dofs(poses);
    CHECK((dofs - returned_dofs).squaredNorm() == Approx(0));
}

// TEST_CASE("Operations on poses", "[physics][pose]")
// {
//     using namespace ccd::physics;
//     int dim = GENERATE(2, 3);
//     int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);
//
//     Eigen::VectorXd sigma_t0 =
//         Eigen::VectorXd::Random(num_bodies * Pose<double>::dim_to_ndof(dim));
//     Eigen::VectorXd sigma_t1 =
//         Eigen::VectorXd::Random(num_bodies * Pose<double>::dim_to_ndof(dim));
//     Eigen::VectorXd delta_sigma = sigma_t1 - sigma_t0;
//
//     Poses<double> poses_t0 = Pose<double>::dofs_to_poses(sigma_t0, dim);
//     Poses<double> poses_t1 = Pose<double>::dofs_to_poses(sigma_t1, dim);
//     Poses<double> delta_poses = poses_t1 - poses_t0;
//     CHECK(
//         Pose<double>::poses_to_dofs(poses_t0 + delta_poses - poses_t1)
//             .squaredNorm()
//         == Approx(0.0).margin(1e-8));
//
//     CHECK(
//         (delta_sigma -
//         Pose<double>::poses_to_dofs(delta_poses)).squaredNorm()
//         == Approx(0.0));
//
//     double t = GENERATE(-1.0, 0.0, 0.5, 0.72, 1.0, 1.24, 3.14);
//     CHECK(
//         (sigma_t0 + delta_sigma * t
//          - Pose<double>::poses_to_dofs(poses_t0 + delta_poses * t))
//             .squaredNorm()
//         == Approx(0.0));
// }

TEST_CASE("Cast poses", "[physics][pose]")
{
    using namespace ccd::physics;
    int dim = GENERATE(2, 3);
    int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);

    Eigen::VectorXf dof =
        Eigen::VectorXf::Random(num_bodies * Pose<float>::dim_to_ndof(dim));
    Poses<float> expected_posesf = Pose<float>::dofs_to_poses(dof, dim);

    Poses<float> actual_posesf =
        cast<double, float>(cast<float, double>(expected_posesf));

    CHECK(
        (Pose<float>::poses_to_dofs(expected_posesf)
         - Pose<float>::poses_to_dofs(actual_posesf))
            .squaredNorm()
        == Approx(0.0));
}
