#include <catch2/catch.hpp>

#include <iostream>
#include <opt/distance_barrier_constraint.hpp>
#include <opt/barrier.hpp>

TEST_CASE("Distance Barrier Constraint", "[opt][ccd][DistanceBarrier]")
{
    using namespace ccd::opt;
    DistanceBarrierConstraint barrier;

    double barrier_epsilon = 1.0;
    barrier.detection_method = ccd::BRUTE_FORCE;
    barrier.custom_inital_epsilon = barrier_epsilon;
    barrier.set_barrier_epsilon(barrier_epsilon);

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(2, 2);
    Eigen::MatrixX2d displacements(4, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    // initial configuration
    //      |3     edge_1
    //      |2
    //
    //   0-------1 edge_0
    //
    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;
    vertices.row(2) << 0.0, 0.5;
    vertices.row(3) << 0.0, 1.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    Eigen::VectorXd expected_barrier = Eigen::VectorXd(8);

    SECTION("No displacements")
    {
        displacements.row(2) << 0.0, 0.0;
        displacements.row(3) << 0.0, 0.0;
        // v0-e0, v0-e1, v1-e0, v1-e1, v2-e0, v2-e1, v3-e0, v3-e1

        // clang-format off
        expected_barrier <<
                0.0, spline_barrier<double>(0.5 * sqrt(2), barrier_epsilon),
                0.0, spline_barrier<double>(0.5 * sqrt(2), barrier_epsilon),
                spline_barrier<double>(0.5, barrier_epsilon), 0.0,
                spline_barrier<double>(1.0, barrier_epsilon), 0.0;
        // clang-format on
    }
    Eigen::VectorXd actual_barrier;
    barrier.initialize(vertices, edges, displacements);
    barrier.compute_constraints(displacements, actual_barrier);
    REQUIRE(actual_barrier.rows() == expected_barrier.rows());
    CHECK((actual_barrier - expected_barrier).squaredNorm() < 1e-16);
}
