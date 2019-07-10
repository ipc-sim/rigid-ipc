#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <autodiff/finitediff.hpp>
#include <opt/rigid_body_problem.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE(
    "Rigid Body Problem Functional", "[RB][RB-Problem][RB-Problem-functional]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d displ_1, displ_2;

    Eigen::MatrixX2d expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    // expected displacement of nodes
    double dx = 0.0;

    SECTION("Translation Case")
    {
        displ_1 << 0.5, 0.5, 0.0;
        displ_2 << 1.0, 1.0, 0.0;
        dx = 0.5 * 0.5 * 8;
    }

    SECTION("90 Deg Rotation Case")
    {
        displ_1 << 0.0, 0.0, 0.5 * M_PI;
        displ_2 << 0.0, 0.0, M_PI;
        dx = 1.0 * 1.0 * 4;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(RigidBody::from_displacement(vertices, edges, displ_1, false));
    rbs.push_back(RigidBody::from_displacement(vertices, edges, displ_2, false));

    auto rbp = RigidBodyProblem();
    auto cstr = BarrierConstraint();
    rbp.init(rbs, cstr);

    // displacement cases
    Eigen::VectorXd x(6);
    x << rbs[0].position_prev + displ_1, rbs[1].position_prev + displ_2;
    double fx = rbp.eval_f(x);
    CHECK(fx == Approx(0.0));

    x << rbs[0].position_prev + 2 * displ_1, rbs[1].position_prev + displ_2;
    fx = rbp.eval_f(x);
    CHECK(fx == Approx(dx / 2));
}

TEST_CASE(
    "Rigid Body Problem Gradient", "[RB][RB-Problem][RB-Problem-gradient]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d vel_1, vel_2;

    Eigen::MatrixX2d expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        vel_1 << 0.5, 0.5, 0.0;
        vel_2 << 1.0, 1.0, 0.0;
    }

    SECTION("90 Deg Rotation Case")
    {
        vel_1 << 0.0, 0.0, 0.5 * M_PI;
        vel_2 << 0.0, 0.0, M_PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1 << 0.5, 0.5, 0.5 * M_PI;
        vel_2 << 1.0, 1.0, M_PI;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(RigidBody::from_displacement(vertices, edges, vel_1, false));
    rbs.push_back(RigidBody::from_displacement(vertices, edges, vel_2, false));

    auto rbp = RigidBodyProblem();
    auto cstr = BarrierConstraint();
    rbp.init(rbs, cstr);


    // displacement cases
    Eigen::VectorXd x(6);
    x << vel_1, vel_2;
    Eigen::VectorXd grad_fx = rbp.eval_grad_f(x);
    Eigen::VectorXd grad_fx_approx = rbp.eval_grad_f_approx(x);

    CHECK(ccd::compare_gradient(grad_fx, grad_fx_approx));
}

TEST_CASE(
    "Rigid Body Problem Hessian", "[RB][RB-Problem][RB-Problem-hessian]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d vel_1, vel_2;

    Eigen::MatrixX2d expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        vel_1 << 0.5, 0.5, 0.0;
        vel_2 << 1.0, 1.0, 0.0;
    }

    SECTION("90 Deg Rotation Case")
    {
        vel_1 << 0.0, 0.0, 0.5 * M_PI;
        vel_2 << 0.0, 0.0, M_PI;
    }
    SECTION("Translation and Rotation Case")
    {
        vel_1 << 0.5, 0.5, 0.5 * M_PI;
        vel_2 << 1.0, 1.0, M_PI;
    }

    using namespace ccd::physics;
    using namespace ccd::opt;

    std::vector<RigidBody> rbs;
    rbs.push_back(RigidBody::from_displacement(vertices, edges, vel_1, false));
    rbs.push_back(RigidBody::from_displacement(vertices, edges, vel_2, false));

    auto rbp = RigidBodyProblem();
    auto cstr = BarrierConstraint();
    rbp.init(rbs, cstr);

    // displacement cases
    Eigen::VectorXd x(6);
    x << vel_1, vel_2;
    Eigen::MatrixXd hess_fx = rbp.eval_hessian_f(x).toDense();
    Eigen::MatrixXd hess_fx_approx = rbp.eval_hessian_f_approx(x);

    CHECK(ccd::compare_jacobian(hess_fx, hess_fx_approx));
}
