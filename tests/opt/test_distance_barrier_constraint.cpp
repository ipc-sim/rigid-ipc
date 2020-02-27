#include <catch2/catch.hpp>

#include <finitediff.hpp>

#include <barrier/barrier.hpp>
#include <logger.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <physics/rigid_body_assembler.hpp>
#include <utils/flatten.hpp>

using namespace ccd;
using namespace opt;
using namespace physics;

TEST_CASE(
    "Distance Barrier Constraint",
    "[opt][ccd][DistanceBarrier][DistanceBarrierConstraint]")
{
    DistanceBarrierConstraint barrier;

    double barrier_epsilon = GENERATE(0.5, 1.0, 5.0);

    barrier.detection_method = ccd::DetectionMethod::BRUTE_FORCE;
    barrier.custom_inital_epsilon = barrier_epsilon;
    barrier.set_barrier_epsilon(barrier_epsilon);
    barrier.active_constraint_scale = 10000;

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Poses<double> displacements(2, Pose<double>::Zero(/*dim=*/2));

    edges.row(0) << 0, 1;
    edges.row(1) << 0, 1;

    // initial configuration
    //      |3     edge_1
    //      |2
    //
    //   0-------1 edge_0
    //
    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;
    vertices.row(2) << 0.0, -0.25;
    vertices.row(3) << 0.0, 0.25;

    Eigen::VectorXd expected_barrier = Eigen::VectorXd(4);
    SECTION("No displacements")
    {
        displacements[1].position << 0.0, 0.0;
        // clang-format off
        // loop over edges then vertices
        expected_barrier <<
            barrier.distance_barrier<double>(0.5, barrier_epsilon),
            barrier.distance_barrier<double>(1.0, barrier_epsilon),
            barrier.distance_barrier<double>(0.5 * sqrt(2),barrier_epsilon),
            barrier.distance_barrier<double>(0.5 * sqrt(2), barrier_epsilon);
        // clang-format on
    }

    SECTION("Left displacements")
    {
        displacements[1].position << -0.5, 0.0;

        // clang-format off
        // loop over edges then vertices
        expected_barrier <<
            barrier.distance_barrier<double>(0.5, barrier_epsilon),
            barrier.distance_barrier<double>(1.0, barrier_epsilon),
            barrier.distance_barrier<double>(0.5, barrier_epsilon),
            barrier.distance_barrier<double>(sqrt(0.5 * 0.5 + 1),
            barrier_epsilon);
        // clang-format on
    }

    SECTION("Farther Left displacements")
    {
        displacements[1].position << -1.0, 0.0;

        // clang-format off
        // loop over edges then vertices
        expected_barrier <<
            barrier.distance_barrier<double>(0.5 * sqrt(2), barrier_epsilon),
            barrier.distance_barrier<double>(sqrt(0.5*0.5 + 1.0),
            barrier_epsilon), barrier.distance_barrier<double>(0.5 * sqrt(2),
            barrier_epsilon), barrier.distance_barrier<double>(sqrt(0.5 * 0.5
            + 1.5 * 1.5), barrier_epsilon);
        // clang-format on
    }

    Eigen::VectorXd actual_barrier;
    // use brute force so we know the order
    RigidBody bodyA = RigidBody::from_points(
        vertices.topRows(2), edges.topRows(1), Pose<double>::Zero(/*dim=*/2),
        displacements[0], 1.0, Eigen::VectorXb::Zero(3), false);
    RigidBody bodyB = RigidBody::from_points(
        vertices.bottomRows(2), edges.bottomRows(1), Pose<double>(0, 0.75, 0),
        displacements[1], 1.0, Eigen::VectorXb::Zero(3), false);

    RigidBodyAssembler rbs;
    rbs.init({ { bodyA, bodyB } });

    barrier.initialize();
    barrier.compute_constraints(
        rbs, rbs.rb_poses(), displacements, actual_barrier);

    REQUIRE(actual_barrier.size() == expected_barrier.size());
    CAPTURE(logger::fmt_eigen(actual_barrier));
    for (int i = 0; i < expected_barrier.size(); i++) {
        CAPTURE(i);
        CHECK(actual_barrier[i] == Approx(expected_barrier[i]));
    }
}

// TEST_CASE(
//     "Distance Barrier Constraint Gradient",
//     "[opt][ccd][DistanceBarrier][DistanceBarrierGradient]")
// {
//     DistanceBarrierConstraint barrier;
//
//     double barrier_epsilon = GENERATE(0.01, 0.5, 1.0, 5.0);
//
//     barrier.detection_method = ccd::DetectionMethod::BRUTE_FORCE;
//     barrier.custom_inital_epsilon = barrier_epsilon;
//     barrier.set_barrier_epsilon(barrier_epsilon);
//     barrier.active_constraint_scale = 10000;
//
//     Eigen::MatrixX2d vertices(4, 2);
//     Eigen::MatrixX2i edges(2, 2);
//     Poses<double> displacements(2, Pose<double>::Zero(/*dim=*/2));
//
//     edges.row(0) << 0, 1;
//     edges.row(1) << 0, 1;
//
//     // initial configuration
//     //      |3     edge_1
//     //      |2
//     //
//     //   0-------1 edge_0
//     //
//     vertices.row(0) << -0.5, 0.0;
//     vertices.row(1) << 0.5, 0.0;
//     vertices.row(2) << 0.0, -0.25;
//     vertices.row(3) << 0.0, 0.25;
//
//     SECTION("No displacements") { displacements[1].position << 0.0, 0.0; }
//
//     SECTION("Left displacements") { displacements[1].position << -0.5, 0.0; }
//
//     SECTION("Farther Left displacements")
//     {
//         displacements[1].position << -1.0, 0.0;
//     }
//
//     RigidBody bodyA = RigidBody::from_points(
//         vertices.topRows(2), edges.topRows(1), Pose<double>::Zero(/*dim=*/2),
//         displacements[0], 1.0, Eigen::VectorXb::Zero(3), false);
//     RigidBody bodyB = RigidBody::from_points(
//         vertices.bottomRows(2), edges.bottomRows(1), Pose<double>(0, 0.75, 0),
//         displacements[1], 1.0, Eigen::VectorXb::Zero(3), false);
//
//     RigidBodyAssembler rbs;
//     rbs.init({ { bodyA, bodyB } });
//
//     Eigen::MatrixXd actual_jac;
//     barrier.initialize();
//     barrier.compute_constraints_jacobian(
//         rbs, rbs.rb_poses(), displacements, actual_jac);
//
//     Eigen::MatrixXd approx_jac;
//     auto f = [&](const Eigen::VectorXd& u) -> Eigen::VectorXd {
//         Eigen::VectorXd fx;
//         Eigen::MatrixXd x = u;
//         ccd::unflatten(x, 2);
//         barrier.compute_constraints(x, fx);
//         barrier.compute_constraints(
//             rbs, rbs.rb_poses(), Pose<double>::dofs_to_poses(sigma, /*dim=*/2),
//             fx);
//         return fx;
//     };
//     Eigen::MatrixXd x = Pose<double>::poses_to_dofs(displacements);
//     fd::finite_jacobian(x, f, approx_jac);
//
//     CAPTURE(approx_jac.rows(), approx_jac.cols());
//     REQUIRE(approx_jac.rows() == actual_jac.rows());
//     REQUIRE(approx_jac.cols() == actual_jac.cols());
//     CHECK((approx_jac - actual_jac).squaredNorm() < 1e-12);
// }
//
// TEST_CASE(
//     "Distance Barrier Constraint Hessian",
//     "[opt][ccd][DistanceBarrier][DistanceBarrierHessian]")
// {
//     DistanceBarrierConstraint barrier;
//
//     double barrier_epsilon = GENERATE(0.5, 1.0, 5.0);
//
//     barrier.detection_method = ccd::DetectionMethod::BRUTE_FORCE;
//     barrier.custom_inital_epsilon = barrier_epsilon;
//     barrier.set_barrier_epsilon(barrier_epsilon);
//
//     Eigen::MatrixX2d vertices(4, 2);
//     Eigen::MatrixX2i edges(2, 2);
//     Poses<double> displacements(2, Pose<double>::Zero(/*dim=*/2));
//
//     edges.row(0) << 0, 1;
//     edges.row(1) << 0, 1;
//
//     // initial configuration
//     //      |3     edge_1
//     //      |2
//     //
//     //   0-------1 edge_0
//     //
//     vertices.row(0) << -0.5, 0.0;
//     vertices.row(1) << 0.5, 0.0;
//     vertices.row(2) << 0.0, -0.25;
//     vertices.row(3) << 0.0, 0.25;
//
//     SECTION("No displacements") { displacements[1].position << 0.0, 0.0; }
//     SECTION("Left displacements") { displacements[1].position << -0.49, 0.0;
//     } SECTION("Farther Left displacements")
//     {
//         displacements[1].position << -1.0, 0.0;
//     }
//     SECTION("Up displacements") { displacements[1].position << 0.0, 0.1; }
//     SECTION("Down displacements") { displacements[1].position << 0.0, -0.1; }
//
//     RigidBody bodyA = RigidBody::from_points(
//         vertices.topRows(2), edges.topRows(1), Pose<double>::Zero(/*dim=*/2),
//         displacements[0], 1.0, Eigen::VectorXb::Zero(3), false);
//     RigidBody bodyB = RigidBody::from_points(
//         vertices.bottomRows(2), edges.bottomRows(1), Pose<double>::Zero(/*dim=*/2),
//         displacements[1], 1.0, Eigen::VectorXb::Zero(3), false);
//     RigidBodyAssembler rbs;
//     rbs.init({ { bodyA, bodyB } });
//     std::vector<Eigen::SparseMatrix<double>> actual_hess;
//     barrier.initialize();
//     barrier.compute_constraints_hessian(
//         rbs, Poses<double>(2, Pose<double>::Zero(/*dim=*/2)), displacements,
//         actual_hess);
//
//     Eigen::MatrixXd actual_jac;
//
//     Eigen::MatrixXd x = Pose<double>::poses_to_dofs(displacements);
//     for (size_t i = 0; i < actual_hess.size(); i++) {
//         Eigen::MatrixXd finite_hess_i;
//         auto f = [&](const Eigen::VectorXd& sigma) -> Eigen::VectorXd {
//             barrier.compute_constraints_jacobian(
//                 rbs, Poses<double>(2, Pose<double>::Zero(/*dim=*/2)),
//                 Pose<double>::dofs_to_poses(sigma, /*dim=*/2), actual_jac);
//             return actual_jac.row(int(i));
//         };
//         fd::finite_jacobian(x, f, finite_hess_i);
//         bool pass =
//             (finite_hess_i - actual_hess[i].toDense()).squaredNorm() < 1e-12;
//
//         CHECK(pass);
//     }
// }
