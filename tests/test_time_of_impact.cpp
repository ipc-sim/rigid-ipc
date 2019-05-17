#include <iostream>

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <ccd/collision_detection.hpp>
#include <ccd/time_of_impact.hpp>

#include "collision_generator.hpp"

using namespace ccd;

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------

/// Compares the time of impact of different implementations
/// against the expected time of impact
/// Functions:
///     - ccd::compute_edge_vertex_time_of_impact(...)
///     - ccd::autodiff::compute_edge_vertex_time_of_impact<double>(...)
///     - ccd::autodiff::compute_edge_vertex_time_of_impact<DScalar>(...)
void check_toi(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
    const double toi_expected);

/// Compares the autodiff gradient vs finite differences
void check_gradient(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk);

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("TimeOfImpact", "[collision_detection][toi][no-grad]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities

    Eigen::Vector2d* V[3] = { &Vi, &Vj, &Vk };
    Eigen::Vector2d* U[3] = { &Ui, &Uj, &Uk };

    for (size_t i = 0; i < 3; i++) {
        V[i]->setZero();
        U[i]->setZero();
    }

    SECTION("PerpendicularImpact") //  (alpha=0.5)
    {
        Vi << -1.0, 0.0;
        Vj << 1.0, 0.0;
        Vk << 0.0, 1.0;

        // touches, intersects, passes-trough
        auto vel = GENERATE(1.0 + 1e-6, 2.0, 4.0);
        // moving direction:
        //  -> ->, -> -, -> <-, - <-, <- <-
        auto j = GENERATE(0, 1, 2, 3, 4);
        // extension, no-deform, compression,
        auto dx = GENERATE(0.5, 0.0, -0.5);

        Uk << 0.0, -(3 - j) * vel / 2.0;
        Ui << -dx, (j - 1.0) * vel / 2.0;
        Uj << dx, (j - 1.0) * vel / 2.0;

        double toi = 1.0 / vel;
        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, toi);

        SECTION("flipped") { check_toi(Vj, Vi, Vk, Uj, Ui, Uk, toi); }
    }

    SECTION("DoubleImpact") // (rotating edge)
    {
        // See fixtures/double-impact.json

        Vi << -1.0, 0.0;
        Vj << 1.0, 0.0;
        Vk << 0.0, 0.5;

        Ui << 1.6730970740318298, 0.8025388419628143;
        Uj << -1.616142749786377, -0.6420311331748962;
        Uk << 0.0, -1.0;

        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, 0.4482900963);

        SECTION("flipped") { check_toi(Vj, Vi, Vk, Uj, Ui, Uk, 0.4482900963); }
    }
}

TEST_CASE("TimeOfImpactGradient", "[collision_detection][toi][gradient]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities

    Eigen::Vector2d* V[3] = { &Vi, &Vj, &Vk };
    Eigen::Vector2d* U[3] = { &Ui, &Uj, &Uk };

    for (size_t i = 0; i < 3; i++) {
        V[i]->setZero();
        U[i]->setZero();
    }

    SECTION("PerpendicularImpact") //  (alpha=0.5)
    {
        Vi << -1.0, 0.0;
        Vj << 1.0, 0.0;
        Vk << 0.0, 1.0;

        SECTION("VertexMoving")
        {
            // touches, intersects, passes-trough
            auto vel = GENERATE(1.0 + 1e-4, 2.0, 4.0);

            Uk << 0.0, -vel;
            Ui << 0.0, 0.0;
            Uj << 0.0, 0.0;

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
        }

        SECTION("BothMoving")
        {
            // touches, intersects, passes-trough
            auto vel = GENERATE(1.0 + 1e-4, 2.0, 4.0);
            // moving direction:
            //  -> ->, -> -, -> <-, - <-, <- <-
            auto j = GENERATE(0, 1, 2, 3, 4);
            Uk << 0.0, -(3 - j) * vel / 2.0;

            Ui << 0.0, (j - 1) * vel / 2.0;
            Uj << 0.0, (j - 1) * vel / 2.0;

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
        }

        SECTION("Extending")
        {
            // edge is compressing / extending
            auto dx = GENERATE(0.5, -0.5);

            Uk << 0.0, -2.0;
            Ui << -dx, 0.0;
            Uj << dx, 0.0;

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
        }

        SECTION("DoubleImpact") // (rotating edge)
        {
            // See fixtures/double-impact.json

            Vi << -1.0, 0.0;
            Vj << 1.0, 0.0;
            Vk << 0.0, 0.5;

            Ui << 1.6730970740318298, 0.8025388419628143;
            Uj << -1.616142749786377, -0.6420311331748962;
            Uk << 0.0, -1.0;

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
        }
    }
}

TEST_CASE("TimeOfImpactRandom", "[collision_detection][toi][random]")
{

    SECTION("TimeOfImpactRigid")
    {
        using namespace ccd::unittests;

        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));
        check_toi(impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj,
            impact.Uk, impact.toi);

        SECTION("flipped")
        {
            check_toi(impact.Vj, impact.Vi, impact.Vk, impact.Uj, impact.Ui,
                impact.Uk, impact.toi);
        }
    }

    // NOTE: we can't tests TOI for non-rigids since the toi we generate
    // might be the second impact.

    SECTION("GradientRigid")
    {
        using namespace ccd::unittests;
        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));
        check_gradient(
            impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj, impact.Uk);
    }

    SECTION("GradientNonRigid")
    {
        using namespace ccd::unittests;
        auto impact = GENERATE(random_impacts(100, /*rigid=*/false));
        check_gradient(
            impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj, impact.Uk);
    }
}

TEST_CASE("TimeOfImpact-TODO", "[collision_detection][toi][!shouldfail]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities

    Eigen::Vector2d* V[3] = { &Vi, &Vj, &Vk };
    Eigen::Vector2d* U[3] = { &Ui, &Uj, &Uk };

    for (size_t i = 0; i < 3; i++) {
        V[i]->setZero();
        U[i]->setZero();
    }

    SECTION("TangentImpact") //  (alpha=0 || alpha = 1)
    {

        Vi << -0.5, 0.0;
        Vj << -1.5, 0.0;
        Vk << 0.5, 0.0;

        // touches, intersects, passes-trough
        auto vel = GENERATE(1.0 + 1e-6, 2.0, 4.0);
        // moving: both (same), ij, both (op), kl, both (same)
        auto j = GENERATE(0, 1, 2, 3, 4);
        // extension, no-deform, compression,
        auto dx = GENERATE(0.5, 0.0, -0.5);

        Uk << -(3 - j) * vel / 2.0, 0.0;

        Ui << (j - 1.0) * vel / 2.0, 0.0;
        // we only move one so we don't change the toi
        Uj << (j - 1.0) * vel / 2.0, dx;

        double toi = 1.0 / vel;

        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, toi);
        SECTION("flipped")
        {
            // change order of edge indices (i.e edge symmetry)
            check_toi(Vj, Vi, Vk, Uj, Ui, Uk, toi);
        }
    }
}

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------
void check_toi(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
    const double toi_expected)
{
    // check original code
    double toi_actual, alpha;
    bool has_collision = ccd::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi_actual, alpha);
    CHECK(has_collision);
    CHECK(toi_expected == Approx(toi_actual));

    // check autodiff code
    toi_actual = -1.0;
    has_collision = ccd::autodiff::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi_actual);
    CHECK(has_collision);
    CHECK(toi_expected == Approx(toi_actual));

    // check with autodiff variables
    DiffScalarBase::setVariableCount(8);

    DVector2 DUi = dvector(0, Ui);
    DVector2 DUj = dvector(2, Uj);
    DVector2 DUk = dvector(4, Uk);
    DScalar(6, 0.0);
    DScalar(7, 0.0);
    DScalar dtoi_actual;

    has_collision = ccd::autodiff::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, DUi, DUj, DUk, dtoi_actual);
    CHECK(has_collision);
    CHECK(toi_expected == Approx(dtoi_actual.getValue()));
}

void check_gradient(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk)
{
    Vector8d grad, grad_fd;
    ccd::autodiff::compute_edge_vertex_time_of_impact_grad(
        Vi, Vj, Vk, Ui, Uj, Uk, grad);
    REQUIRE(!std::isnan(grad.squaredNorm()));

    ccd::autodiff::compute_edge_vertex_time_of_impact_grad_fd(
        Vi, Vj, Vk, Ui, Uj, Uk, grad_fd);

    REQUIRE(!std::isnan(grad_fd.squaredNorm()));

    bool compare_grad = compare_gradient(grad, grad_fd);
    if (!compare_grad) {
        std::cout << " grad    " << grad.transpose() << std::endl;
        std::cout << " grad_fd " << grad_fd.transpose() << std::endl;
    }
    CHECK(compare_grad);
}
