#include <catch2/catch.hpp>

#include <autodiff/finitediff.hpp>
#include <iostream>
#include <opt/barrier.hpp>
#include <opt/distance_barrier_constraint.hpp>
#include <utils/flatten.hpp>

TEST_CASE("Distance Barrier Constraint", "[opt][ccd][DistanceBarrier]")
{
    using namespace ccd::opt;
    DistanceBarrierConstraint barrier;

    double barrier_epsilon = GENERATE(0.01, 0.5, 1.0, 5.0);

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
                    0.0, spline_barrier<double>(0.5 * sqrt(2),
                    barrier_epsilon), 0.0, spline_barrier<double>(0.5 *
                    sqrt(2), barrier_epsilon), spline_barrier<double>(0.5,
                    barrier_epsilon), 0.0, spline_barrier<double>(1.0,
                    barrier_epsilon), 0.0;
        // clang-format on
    }

    SECTION("Left displacements")
    {
        displacements.row(2) << -0.5, 0.0;
        displacements.row(3) << -0.5, 0.0;
        // v0-e0, v0-e1, v1-e0, v1-e1, v2-e0, v2-e1, v3-e0, v3-e1

        // clang-format off
            expected_barrier <<
                    0.0, spline_barrier<double>(0.5, barrier_epsilon),
                    0.0, spline_barrier<double>(sqrt(0.5 * 0.5 + 1),
                    barrier_epsilon), spline_barrier<double>(0.5,
                    barrier_epsilon), 0.0, spline_barrier<double>(1.0,
                    barrier_epsilon), 0.0;
        // clang-format on
    }

    SECTION("Farther Left displacements")
    {
        displacements.row(2) << -1.0, 0.0;
        displacements.row(3) << -1.0, 0.0;
        // v0-e0, v0-e1, v1-e0, v1-e1, v2-e0, v2-e1, v3-e0, v3-e1

        // clang-format off
            expected_barrier <<
                    0.0, spline_barrier<double>(0.5 * sqrt(2),
                    barrier_epsilon), 0.0, spline_barrier<double>(sqrt(0.5 *
                    0.5 + 1.5 * 1.5), barrier_epsilon),
                    spline_barrier<double>(0.5 * sqrt(2), barrier_epsilon),
                    0.0, spline_barrier<double>(sqrt(0.5*0.5 + 1.0),
                    barrier_epsilon), 0.0;
        // clang-format on
    }
    SECTION("Down displacements, Touchin")
    {
        displacements.row(2) << 0.0, -0.5;
        displacements.row(3) << 0.0, -0.5;
        // v0-e0, v0-e1, v1-e0, v1-e1, v2-e0, v2-e1, v3-e0, v3-e1

        // clang-format off
        expected_barrier <<
                0.0, spline_barrier<double>(0.5, barrier_epsilon),
                0.0, spline_barrier<double>(0.5, barrier_epsilon),
                std::numeric_limits<double>::infinity(), 0.0,
                spline_barrier<double>(0.5, barrier_epsilon), 0.0;
        // clang-format on
    }

    Eigen::VectorXd actual_barrier;
    barrier.initialize(vertices, edges, displacements);
    barrier.compute_constraints(displacements, actual_barrier);
    REQUIRE(actual_barrier.rows() == expected_barrier.rows());
    for (int i = 0; i < expected_barrier.rows(); i++) {
        if (std::isinf(expected_barrier[i]) || std::isinf(actual_barrier[i])) {
            bool match_inf = std::isinf(expected_barrier[i])
                && std::isinf(actual_barrier[i]);
            CHECK(match_inf);
        } else {
            CHECK(actual_barrier[i] == Approx(expected_barrier[i]));
        }
    }
}

TEST_CASE("Distance Barrier Constraint Gradient",
    "[opt][ccd][DistanceBarrier][DistanceBarrierGradient]")
{
    using namespace ccd::opt;
    DistanceBarrierConstraint barrier;

    double barrier_epsilon = GENERATE(0.01, 0.5, 1.0, 5.0);

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

    SECTION("No displacements")
    {
        displacements.row(2) << 0.0, 0.0;
        displacements.row(3) << 0.0, 0.0;
    }

    SECTION("Left displacements")
    {
        displacements.row(2) << -0.5, 0.0;
        displacements.row(3) << -0.5, 0.0;
    }

    SECTION("Farther Left displacements")
    {
        displacements.row(2) << -1.0, 0.0;
        displacements.row(3) << -1.0, 0.0;
    }

    Eigen::MatrixXd actual_jac;
    barrier.initialize(vertices, edges, displacements);
    barrier.compute_constraints_jacobian(displacements, actual_jac);

    Eigen::MatrixXd approx_jac;
    auto f = [&](const Eigen::VectorXd& u) -> Eigen::VectorXd {
        Eigen::VectorXd fx;
        Eigen::MatrixXd x = u;
        ccd::unflatten(x, 2);
        barrier.initialize(vertices, edges, x);
        barrier.compute_constraints(x, fx);
        return fx;
    };
    Eigen::MatrixXd x = displacements;
    ccd::flatten(x);
    ccd::finite_jacobian(x, f, approx_jac);

    CHECK((approx_jac - actual_jac).squaredNorm() < 1e-12);
}

TEST_CASE("Distance Barrier Constraint Hessian",
    "[opt][ccd][DistanceBarrier][DistanceBarrierHessian]")
{
    using namespace ccd::opt;
    DistanceBarrierConstraint barrier;

    double barrier_epsilon = GENERATE(0.01, 0.5, 1.0, 5.0);

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

    SECTION("No displacements")
    {
        displacements.row(2) << 0.0, 0.0;
        displacements.row(3) << 0.0, 0.0;
    }

    SECTION("Left displacements")
    {
        displacements.row(2) << -0.5, 0.0;
        displacements.row(3) << -0.5, 0.0;
    }

    SECTION("Farther Left displacements")
    {
        displacements.row(2) << -1.0, 0.0;
        displacements.row(3) << -1.0, 0.0;
    }

    std::vector<Eigen::SparseMatrix<double>> actual_hess;
    barrier.initialize(vertices, edges, displacements);
    barrier.compute_constraints_hessian(displacements, actual_hess);

    Eigen::MatrixXd actual_jac;

    Eigen::MatrixXd x = displacements;
    ccd::flatten(x);
    for (size_t i = 0; i < actual_hess.size(); i++) {
        Eigen::MatrixXd finite_hess_i;
        auto f = [&](const Eigen::VectorXd& u) -> Eigen::VectorXd {
            Eigen::MatrixXd x = u;
            ccd::unflatten(x, 2);
            barrier.initialize(vertices, edges, x);
            barrier.compute_constraints_jacobian(displacements, actual_jac);
            return actual_jac.row(int(i));
        };

        ccd::finite_jacobian(x, f, finite_hess_i);
        CHECK((finite_hess_i - actual_hess[i].toDense()).squaredNorm() < 1e-12);
    }
}
