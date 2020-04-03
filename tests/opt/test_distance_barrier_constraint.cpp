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
    Poses<double> poses_t1(2, Pose<double>::Zero(/*dim=*/2));

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
        poses_t1[1].position << 0.0, 0.75;
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
        poses_t1[1].position << -0.5, 0.75;

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
        poses_t1[1].position << -1.0, 0.75;

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
        Pose<double>::Zero(/*dim=*/2), 1.0, Eigen::VectorXb::Zero(3), false);
    RigidBody bodyB = RigidBody::from_points(
        vertices.bottomRows(2), edges.bottomRows(1), Pose<double>(0, 0.75, 0),
        Pose<double>::Zero(/*dim=*/2), 1.0, Eigen::VectorXb::Zero(3), false);

    RigidBodyAssembler rbs;
    rbs.init({ { bodyA, bodyB } });

    barrier.initialize();
    barrier.compute_constraints(rbs, poses_t1, actual_barrier);

    REQUIRE(actual_barrier.size() == expected_barrier.size());
    CAPTURE(logger::fmt_eigen(actual_barrier));
    for (int i = 0; i < expected_barrier.size(); i++) {
        CAPTURE(i);
        CHECK(actual_barrier[i] == Approx(expected_barrier[i]));
    }
}
