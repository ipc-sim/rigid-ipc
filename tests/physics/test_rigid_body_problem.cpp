#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <finitediff.hpp>
#include <problems/distance_barrier_rb_problem.hpp>
#include <utils/not_implemented_error.hpp>

namespace test_utils {

using namespace ccd::physics;
RigidBody rb_from_displacements(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    const Pose<double>& displacement)
{
    // move vertices so they center of mass is at 0,0
    Eigen::RowVectorXd x = compute_center_of_mass(vertices, edges);
    Eigen::MatrixXd centered_vertices = vertices.rowwise() - x;

    // set position so current vertices match input
    Pose<double> pose_t0 = Pose<double>::Zero(vertices.cols());
    pose_t0.position = x;

    // set previous_step position to:
    Pose<double> pose_t1 = pose_t0 + displacement;

    // set velocity to zero
    Pose<double> velocity = Pose<double>::Zero(vertices.cols());

    int ndof = pose_t0.ndof();
    auto rb = RigidBody::from_points(
        vertices, edges, pose_t0, velocity,
        /*density=*/1,
        /*dof=*/Eigen::VectorXb::Zero(ndof),
        /*oriented=*/false);
    rb.pose = pose_t1;
    return rb;
}

} // namespace test_utils

TEST_CASE(
    "Rigid Body Problem Functional", "[RB][RB-Problem][RB-Problem-functional]")
{
    using namespace test_utils;

    Eigen::MatrixXd vertices(4, 2);
    int dim = vertices.cols();
    Eigen::MatrixXi edges(4, 2);
    Pose<double> displ_1 = Pose<double>::Zero(dim),
                 displ_2 = Pose<double>::Zero(dim);

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    // expected displacement of nodes
    double dx = 0.0;

    double mass = 4.0; // ∑ mᵢ = 4 * 1.0
    double moment_inertia =
        8.0 * 0.5 * 0.5; // ∑ mᵢ ||rᵢ||² =  4 * 1.0 * 2.0 * (0.5)**2

    SECTION("Translation Case")
    {
        displ_1.position << 0.5, 0.5;
        displ_2.position << 1.0, 1.0;
        dx = 0.5 * mass
            * (displ_1.position.squaredNorm() + displ_2.position.squaredNorm());
    }

    SECTION("90 Deg Rotation Case")
    {
        displ_1.rotation << 0.5 * M_PI;
        displ_2.rotation << M_PI;
        dx = 0.5 * moment_inertia
            * (displ_1.rotation.squaredNorm() + displ_2.rotation.squaredNorm());
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(rb_from_displacements(vertices, edges, displ_1));
    rbs.push_back(rb_from_displacements(vertices, edges, displ_2));

    DistanceBarrierRBProblem rbp("rb_problem");
    rbp.init(rbs);

    // displacement cases
    Eigen::VectorXd x = rbp.poses_to_dofs<double>(
        { { rbs[0].pose_prev + displ_1, rbs[1].pose_prev + displ_2 } });
    double fx = rbp.eval_f(x);
    CHECK(fx == Approx(0.0));

    x = rbp.poses_to_dofs<double>(
        { { rbs[0].pose_prev + displ_1 * 2, rbs[1].pose_prev + displ_2 * 2 } });
    fx = rbp.eval_f(x);
    CHECK(fx == Approx(dx));
}

TEST_CASE(
    "Rigid Body Problem Gradient", "[RB][RB-Problem][RB-Problem-gradient]")
{
    using namespace test_utils;

    Eigen::MatrixXd vertices(4, 2);
    int dim = vertices.cols();
    Eigen::MatrixXi edges(4, 2);
    Pose<double> vel_1 = Pose<double>::Zero(dim),
                 vel_2 = Pose<double>::Zero(dim);

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_2.position << 1.0, 1.0;
    }

    SECTION("90 Deg Rotation Case")
    {
        vel_1.rotation << 0.5 * M_PI;
        vel_2.rotation << M_PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * M_PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << M_PI;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(rb_from_displacements(vertices, edges, vel_1));
    rbs.push_back(rb_from_displacements(vertices, edges, vel_2));

    DistanceBarrierRBProblem rbp("rb_problem");
    rbp.init(rbs);

    // displacement cases
    Eigen::VectorXd x(2 * Pose<double>::dim_to_ndof(dim));
    x << vel_1.dof(), vel_2.dof();
    Eigen::VectorXd grad_fx = rbp.eval_grad_f(x);
    Eigen::VectorXd grad_fx_approx = eval_grad_f_approx(rbp, x);

    CHECK(fd::compare_gradient(grad_fx, grad_fx_approx));
}

TEST_CASE("Rigid Body Problem Hessian", "[RB][RB-Problem][RB-Problem-hessian]")
{
    using namespace test_utils;
    Eigen::MatrixXd vertices(4, 2);
    int dim = vertices.cols();
    Eigen::MatrixXi edges(4, 2);
    Pose<double> vel_1 = Pose<double>::Zero(dim),
                 vel_2 = Pose<double>::Zero(dim);

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_2.position << 1.0, 1.0;
    }

    SECTION("90 Deg Rotation Case")
    {
        vel_1.rotation << 0.5 * M_PI;
        vel_2.rotation << M_PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * M_PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << M_PI;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(rb_from_displacements(vertices, edges, vel_1));
    rbs.push_back(rb_from_displacements(vertices, edges, vel_2));

    DistanceBarrierRBProblem rbp("rb_problem");
    rbp.init(rbs);

    // displacement cases
    Eigen::VectorXd x(2 * Pose<double>::dim_to_ndof(dim));
    x << vel_1.dof(), vel_2.dof();
    Eigen::MatrixXd hess_fx = rbp.eval_hessian_f(x).toDense();
    Eigen::MatrixXd hess_fx_approx = eval_hess_f_approx(rbp, x);

    CHECK(fd::compare_jacobian(hess_fx, hess_fx_approx));
}

TEST_CASE("dof -> poses -> dof", "[RB][RB-Problem]")
{
    using namespace test_utils;
    Eigen::MatrixXd vertices(4, 2);
    int dim = vertices.cols();
    Eigen::MatrixXi edges(4, 2);
    Pose<double> vel_1 = Pose<double>::Zero(dim),
                 vel_2 = Pose<double>::Zero(dim);

    Eigen::MatrixXd expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_2.position << 1.0, 1.0;
    }

    SECTION("90 Deg Rotation Case")
    {
        vel_1.rotation << 0.5 * M_PI;
        vel_2.rotation << M_PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * M_PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << M_PI;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs = {
        { rb_from_displacements(vertices, edges, vel_1),
          rb_from_displacements(vertices, edges, vel_2) }
    };

    DistanceBarrierRBProblem rbp("rb_problem");
    rbp.init(rbs);

    for (int i = 0; i < 100; i++) {
        Eigen::VectorXd expected_dof = Eigen::VectorXd::Random(3 * rbs.size());
        Eigen::VectorXd actual_dof =
            rbp.poses_to_dofs(rbp.dofs_to_poses(expected_dof));
        REQUIRE(expected_dof.size() == actual_dof.size());
        CHECK((expected_dof - actual_dof).squaredNorm() == Approx(0.0));
    }
}

// TODO: Add 3D RB test
