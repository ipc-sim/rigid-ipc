#include <iostream>

#include <catch2/catch.hpp>

#include <ccd/ccd.hpp>
#include <ccd/linear/edge_vertex_ccd.hpp>
#include <finitediff.hpp>

#include "collision_generator.hpp"

using namespace ipc::rigid;

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------

/// Compares the time of impact of different implementations
/// against the expected time of impact
/// Functions:
///     - ipc::rigid::compute_edge_vertex_time_of_impact(...)
void check_toi(
    const Eigen::Vector2d& Vi,
    const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk,
    const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj,
    const Eigen::Vector2d& Uk,
    const double toi_expected)
{
    double toi_actual;
    // check autodiff code
    toi_actual = -1.0;
    bool has_collision = ipc::rigid::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi_actual);
    CHECK(has_collision);
    CHECK(toi_expected == Approx(toi_actual));
}

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

TEST_CASE("TimeOfImpactRandom", "[collision_detection][toi][random]")
{

    SECTION("TimeOfImpactRigid")
    {
        using namespace ipc::rigid::unittests;

        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));
        check_toi(
            impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj, impact.Uk,
            impact.toi);

        SECTION("flipped")
        {
            check_toi(
                impact.Vj, impact.Vi, impact.Vk, impact.Uj, impact.Ui,
                impact.Uk, impact.toi);
        }
    }
}

TEST_CASE("TimeOfImpactBadCases", "[collision_detection][toi]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities

    Vi = Vj = Vk = Ui = Uj = Uk = Eigen::Vector2d::Zero();

    SECTION("TangentImpact") //  (alpha=0 || alpha = 1)
    {
        Vi << -0.5, 0.0;
        Vj << -1.5, 0.0;
        Vk << 0.5, 0.0;

        // touches, intersects, passes-trough
        auto vel = GENERATE(1.0 + 1e-3, 2.0, 4.0);
        // moving: both (same), ij, both (op), kl, both (same)
        auto j = GENERATE(0, 1, 2, 3, 4);
        // extension, no-deform, compression,
        auto dx = GENERATE(0.5, 0.0, -0.5);

        Uk << -(3 - j) * vel / 2.0, 0.0;

        Ui << (j - 1.0) * vel / 2.0, 0.0;
        // we only move one so we don't change the toi
        Uj << (j - 1.0) * vel / 2.0, dx;

        double toi = 1.0 / vel;

        CAPTURE(vel, j, dx);
        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, toi);
        // change order of edge indices (i.e edge symmetry)
        check_toi(Vj, Vi, Vk, Uj, Ui, Uk, toi);
    }
}
