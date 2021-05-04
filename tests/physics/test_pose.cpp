// Test Pose
#include <catch2/catch.hpp>

#include <Eigen/Geometry>
#include <igl/PI.h>

#include <autodiff/autodiff_types.hpp>
#include <interval/interval.hpp>
#include <physics/pose.hpp>

TEST_CASE("Poses to dofs", "[physics][pose]")
{
    using namespace ipc::rigid;
    int dim = GENERATE(2, 3);
    int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);
    Eigen::VectorXd dofs =
        Eigen::VectorXd::Random(num_bodies * Pose<double>::dim_to_ndof(dim));
    Poses<double> poses = Pose<double>::dofs_to_poses(dofs, dim);
    Eigen::VectorXd returned_dofs = Pose<double>::poses_to_dofs(poses);
    CHECK((dofs - returned_dofs).squaredNorm() == Approx(0));
}

TEST_CASE("Cast poses", "[physics][pose]")
{
    using namespace ipc::rigid;
    int dim = GENERATE(2, 3);
    int num_bodies = GENERATE(0, 1, 2, 3, 10, 1000);

    Eigen::VectorXf dof =
        Eigen::VectorXf::Random(num_bodies * Pose<float>::dim_to_ndof(dim));
    Poses<float> expected_posesf = Pose<float>::dofs_to_poses(dof, dim);

    Poses<float> actual_posesf = cast<float>(cast<double>(expected_posesf));

    CHECK(
        (Pose<float>::poses_to_dofs(expected_posesf)
         - Pose<float>::poses_to_dofs(actual_posesf))
            .squaredNorm()
        == Approx(0.0));
}

TEST_CASE("SE(3) ↦ SO(3)", "[physics][pose]")
{
    using namespace ipc::rigid;
    double angle;
    Eigen::Vector3d axis;

    SECTION("zero")
    {
        angle = 0;
        axis = Eigen::Vector3d::Random();
    }
    SECTION("random")
    {
        angle = GENERATE(take(100, random(0.0, 2 * igl::PI)));
        axis = Eigen::Vector3d::Random();
    }
    axis.normalize();

    Pose<double> p = Pose<double>::Zero(3);
    p.rotation = angle * axis;
    Eigen::Matrix3d R_actual = p.construct_rotation_matrix();
    Eigen::Matrix3d R_expected =
        Eigen::AngleAxisd(angle, axis).toRotationMatrix();
    CHECK((R_actual - R_expected).norm() == Approx(0).margin(1e-12));
}

TEST_CASE("∇²(SE(3) ↦ SO(3))", "[!benchmark][physics][pose]")
{
    using namespace ipc::rigid;
    typedef ipc::rigid::AutodiffType<Eigen::Dynamic, 12> Diff;
    Diff::activate(12);

    Pose<double> p;
    p.position = Eigen::Vector3d::Zero();
    p.rotation = Eigen::Vector3d(0, igl::PI, 0);

    BENCHMARK("Compute R") { return p.construct_rotation_matrix(); };

    Pose<Diff::DDouble1> d1p;
    d1p.position = Diff::d1vars(0, Eigen::Vector3d::Zero());
    d1p.rotation = Diff::d1vars(3, Eigen::Vector3d(0, igl::PI, 0));

    BENCHMARK("Compute R DDouble1") { return d1p.construct_rotation_matrix(); };

    Pose<Diff::DDouble2> d2p;
    d2p.position = Diff::d2vars(0, Eigen::Vector3d::Zero());
    d2p.rotation = Diff::d2vars(3, Eigen::Vector3d(0, igl::PI, 0));

    BENCHMARK("Compute R DDouble2") { return d2p.construct_rotation_matrix(); };
}

TEST_CASE("Interval SE(3) ↦ SO(3)", "[!benchmark][physics][pose]")
{
    using namespace ipc::rigid;
    double angle;
    Eigen::Vector3d axis;

    SECTION("zero")
    {
        angle = 0;
        axis = Eigen::Vector3d::Random();
    }
    SECTION("random")
    {
        angle = GENERATE(take(1, random(0.0, 2 * igl::PI)));
        axis = Eigen::Vector3d::Random();
    }
    axis.normalize();

    Pose<double> p = Pose<double>::Zero(3);
    p.rotation = angle * axis;
    BENCHMARK("Double SE(3) ↦ SO(3)")
    {
        Eigen::Matrix3d R = p.construct_rotation_matrix();
    };
    Pose<ipc::rigid::Interval> pI = p.cast<ipc::rigid::Interval>();
    BENCHMARK("Interval SE(3) ↦ SO(3)")
    {
        Matrix3I R = pI.construct_rotation_matrix();
    };
}
