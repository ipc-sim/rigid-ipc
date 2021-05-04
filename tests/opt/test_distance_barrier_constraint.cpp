#include <catch2/catch.hpp>

#include <logger.hpp>
#include <opt/distance_barrier_constraint.hpp>

using namespace ipc;
using namespace ipc::rigid;

// TEST_CASE(
//     "Distance Barrier Constraint",
//     "[opt][ccd][DistanceBarrier][DistanceBarrierConstraint]")
// {
//     DistanceBarrierConstraint barrier;
//
//     double barrier_epsilon = GENERATE(0.5, 1.0, 5.0);
//
//     barrier.detection_method = ipc::rigid::DetectionMethod::BRUTE_FORCE;
//     barrier.initial_barrier_activation_distance = barrier_epsilon;
//     barrier.barrier_activation_distance(barrier_epsilon);
//     barrier.active_constraint_scale = 10000;
//
//     Eigen::MatrixX2d vertices(4, 2);
//     Eigen::MatrixX2i edges(2, 2);
//     Poses<double> poses_t1(2, Pose<double>::Zero(/*dim=*/2));
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
//     auto expected_barrier_func = [&](double d) {
// #ifdef USE_DISTANCE_SQUARED
//         return barrier.distance_barrier<double>(d * d, barrier_epsilon);
// #else
//         return barrier.distance_barrier<double>(d, barrier_epsilon);
// #endif
//     };
//
//     Eigen::VectorXd expected_barrier = Eigen::VectorXd(4);
//     SECTION("No displacements")
//     {
//         poses_t1[1].position << 0.0, 0.75;
//         expected_barrier << expected_barrier_func(0.5),
//             expected_barrier_func(1.0), expected_barrier_func(0.5 * sqrt(2)),
//             expected_barrier_func(0.5 * sqrt(2));
//     }
//
//     SECTION("Left displacements")
//     {
//         poses_t1[1].position << -0.5, 0.75;
//
//         expected_barrier << expected_barrier_func(0.5),
//             expected_barrier_func(1.0), expected_barrier_func(0.5),
//             expected_barrier_func(sqrt(0.5 * 0.5 + 1));
//     }
//
//     SECTION("Farther Left displacements")
//     {
//         poses_t1[1].position << -1.0, 0.75;
//
//         expected_barrier << expected_barrier_func(0.5 * sqrt(2)),
//             expected_barrier_func(sqrt(0.5 * 0.5 + 1.0)),
//             expected_barrier_func(0.5 * sqrt(2)),
//             expected_barrier_func(sqrt(0.5 * 0.5 + 1.5 * 1.5));
//     }
//
//     Eigen::VectorXd actual_barrier;
//     // use brute force so we know the order
//     RigidBody bodyA = RigidBody(
//         vertices.topRows(2), edges.topRows(1), Pose<double>::Zero(/*dim=*/2),
//         /*velocity=*/Pose<double>::Zero(/*dim=*/2),
//         /*force=*/Pose<double>::Zero(/*dim=*/2), /*density=*/1.0,
//         /*is_dof_fixed=*/VectorXb::Zero(3), /*is_oriented=*/false,
//         /*group=*/0);
//     RigidBody bodyB = RigidBody(
//         vertices.bottomRows(2), edges.bottomRows(1), Pose<double>(0, 0.75,
//         0),
//         /*velocity=*/Pose<double>::Zero(/*dim=*/2),
//         /*force=*/Pose<double>::Zero(/*dim=*/2), /*density=*/1.0,
//         /*is_dof_fixed=*/VectorXb::Zero(3), /*is_oriented=*/false,
//         /*group=*/1);
//
//     RigidBodyAssembler rbs;
//     rbs.init({ { bodyA, bodyB } });
//
//     barrier.initialize();
//     barrier.compute_constraints(rbs, poses_t1, actual_barrier);
//
//     REQUIRE(actual_barrier.size() == expected_barrier.size());
//     CAPTURE(logger::fmt_eigen(actual_barrier));
//     for (int i = 0; i < expected_barrier.size(); i++) {
//         CAPTURE(i);
//         CHECK(actual_barrier[i] == Approx(expected_barrier[i]));
//     }
// }
