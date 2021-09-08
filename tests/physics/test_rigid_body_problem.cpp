#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>
#include <finitediff.hpp>
#include <igl/PI.h>

#include <physics/mass.hpp>
#include <problems/split_distance_barrier_rb_problem.hpp>
#include <utils/not_implemented_error.hpp>

using namespace ipc;
using namespace ipc::rigid;

RigidBody rb_from_displacements(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    Pose<double> pose_t1)
{
    static int id = 0;
    // move vertices so they center of mass is at 0,0
    VectorMax3d x = compute_center_of_mass(vertices, edges);
    Eigen::MatrixXd centered_vertices = vertices.rowwise() - x.transpose();

    // set position so current vertices match input
    Pose<double> pose_t0 = Pose<double>::Zero(vertices.cols());
    pose_t0.position = x;

    // set previous_step position to:
    pose_t1.position += pose_t0.position;

    int ndof = pose_t0.ndof();
    auto rb = RigidBody(
        vertices, edges, pose_t0,
        /*velocity=*/Pose<double>::Zero(vertices.cols()),
        /*force=*/Pose<double>::Zero(vertices.cols()),
        /*density=*/1,
        /*dof=*/VectorXb::Zero(ndof),
        /*oriented=*/false,
        /*group=*/id++);
    rb.pose = pose_t1;
    return rb;
}

TEST_CASE(
    "Rigid Body Problem Functional", "[RB][RB-Problem][RB-Problem-functional]")
{
    Eigen::MatrixXd vertices(4, 2);
    int dim = vertices.cols();
    Eigen::MatrixXi edges(4, 2);
    Pose<double> rb1_pose_t1 = Pose<double>::Zero(dim);
    Pose<double> rb2_pose_t1 = Pose<double>::Zero(dim);

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
        rb1_pose_t1.position << 0.5, 0.5;
        rb2_pose_t1.position << 1.0, 1.0;
        dx = 0.5 * mass
            * (rb1_pose_t1.position.squaredNorm()
               + rb2_pose_t1.position.squaredNorm());
    }

    SECTION("90 Deg Rotation Case")
    {
        rb1_pose_t1.rotation << 0.5 * igl::PI;
        rb2_pose_t1.rotation << igl::PI;
        dx = 0.5 * moment_inertia
            * (rb1_pose_t1.rotation.squaredNorm()
               + rb2_pose_t1.rotation.squaredNorm());
    }

    std::vector<RigidBody> rbs = {
        { rb_from_displacements(vertices, edges, rb1_pose_t1),
          rb_from_displacements(vertices, edges, rb2_pose_t1) }
    };

    SplitDistanceBarrierRBProblem rbp;
    rbp.init(rbs);

    // displacement cases

    Eigen::VectorXd x =
        rbp.poses_to_dofs<double>({ { rb1_pose_t1, rb2_pose_t1 } });
    double fx = rbp.compute_energy_term(x);
    CHECK(fx == Approx(0.0));

    x = rbp.poses_to_dofs<double>({ { 2 * rb1_pose_t1, 2 * rb2_pose_t1 } });
    fx = rbp.compute_energy_term(x);
    CHECK(fx == Approx(dx));
}

TEST_CASE(
    "Rigid Body Problem Gradient", "[RB][RB-Problem][RB-Problem-gradient]")
{
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
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.rotation << igl::PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << igl::PI;
    }

    std::vector<RigidBody> rbs;
    rbs.push_back(rb_from_displacements(vertices, edges, vel_1));
    rbs.push_back(rb_from_displacements(vertices, edges, vel_2));

    SplitDistanceBarrierRBProblem rbp;
    rbp.init(rbs);

    // displacement cases
    Eigen::VectorXd x(2 * Pose<double>::dim_to_ndof(dim));
    x << vel_1.dof(), vel_2.dof();
    Eigen::VectorXd grad_fx;
    double fx = rbp.compute_energy_term(x, grad_fx);
    Eigen::VectorXd grad_fx_approx = eval_grad_energy_approx(rbp, x);

    CHECK(fd::compare_gradient(grad_fx, grad_fx_approx));
}

TEST_CASE("Rigid Body Problem Hessian", "[RB][RB-Problem][RB-Problem-hessian]")
{
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
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.rotation << igl::PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << igl::PI;
    }

    std::vector<RigidBody> rbs;
    rbs.push_back(rb_from_displacements(vertices, edges, vel_1));
    rbs.push_back(rb_from_displacements(vertices, edges, vel_2));

    SplitDistanceBarrierRBProblem rbp;
    rbp.init(rbs);

    // displacement cases
    Eigen::VectorXd x(2 * Pose<double>::dim_to_ndof(dim));
    x << vel_1.dof(), vel_2.dof();

    Eigen::VectorXd grad_fx;
    Eigen::SparseMatrix<double> hess_fx;
    rbp.compute_energy_term(x, grad_fx, hess_fx);
    Eigen::MatrixXd hess_fx_approx = eval_hess_energy_approx(rbp, x);

    CHECK(fd::compare_jacobian(hess_fx.toDense(), hess_fx_approx));
}

TEST_CASE("dof -> poses -> dof", "[RB][RB-Problem]")
{
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
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.rotation << igl::PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1.position << 0.5, 0.5;
        vel_1.rotation << 0.5 * igl::PI;
        vel_2.position << 1.0, 1.0;
        vel_2.rotation << igl::PI;
    }

    std::vector<RigidBody> rbs = {
        { rb_from_displacements(vertices, edges, vel_1),
          rb_from_displacements(vertices, edges, vel_2) }
    };

    SplitDistanceBarrierRBProblem rbp;
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
