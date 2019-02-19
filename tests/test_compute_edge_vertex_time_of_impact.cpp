#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch.hpp>

#include <FixingCollisions/collision_detection.hpp>
#include <FixingCollisions/degenerate_edge_error.hpp>
#include <autodiff/finitediff.hpp>
#include <autogen/collision_volume.hpp>

using namespace ccd;

void check_toi(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
    const double toi_expected, const bool v2 = false)
{
    double toi_actual, alpha;
    bool has_collision = ccd::compute_edge_vertex_time_of_impact(
        Vk, Uk, Vi, Ui, Vj, Uj, toi_actual, alpha);
    REQUIRE(has_collision);
    REQUIRE(toi_expected == Approx(toi_actual));

    toi_actual = -1.0;
    has_collision = ccd::autogen::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, Ui, Uj, Uk, toi_actual);
    REQUIRE(has_collision);
    REQUIRE(toi_expected == Approx(toi_actual));

    if (v2) {
        toi_actual = -1.0;
        has_collision = ccd::autogen::compute_edge_vertex_time_of_impact_v2(
            Vi, Vj, Vk, Ui, Uj, Uk, toi_actual);

        REQUIRE(has_collision);
        REQUIRE(toi_expected == Approx(toi_actual));
    }

    // check with autodiff variables
    DiffScalarBase::setVariableCount(8);

    DVector2 DUi = dvector(0, Ui);
    DVector2 DUj = dvector(2, Uj);
    DVector2 DUk = dvector(4, Uk);
    DScalar(6, 0.0);
    DScalar(7, 0.0);
    DScalar dtoi_actual;

    has_collision = ccd::autogen::compute_edge_vertex_time_of_impact(
        Vi, Vj, Vk, DUi, DUj, DUk, dtoi_actual);
    REQUIRE(has_collision);
    REQUIRE(toi_expected == Approx(dtoi_actual.getValue()));
}

void check_gradient(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk)
{
    Vector8d grad, grad_fd;
    ccd::autogen::compute_edge_vertex_time_of_impact_grad(
        Vi, Vj, Vk, Ui, Uj, Uk, /*version=*/2, grad);
    REQUIRE(!std::isnan(grad.squaredNorm()));

    ccd::autogen::compute_edge_vertex_time_of_impact_grad_fd(
        Vi, Vj, Vk, Ui, Uj, Uk, grad_fd);

    REQUIRE(!std::isnan(grad_fd.squaredNorm()));

    REQUIRE(compare_gradient(grad, grad_fd));
}
TEST_CASE("TestTOI", "[toi][collision_detection]")
{
    Eigen::Vector2d Vi, Vj, Vk, Vl; // positions
    Eigen::Vector2d Ui, Uj, Uk, Ul; // velocities

    Eigen::Vector2d* V[4] = { &Vi, &Vj, &Vk, &Vl };
    Eigen::Vector2d* U[4] = { &Ui, &Uj, &Uk, &Ul };

    for (size_t i = 0; i < 4; i++) {
        V[i]->setZero();
        U[i]->setZero();
    }

    SECTION("PerpendicularImpact") //  (alpha=0.5)
    {
        Vi << -1.0, 0.0;
        Vj << 1.0, 0.0;
        Vk << 0.0, 1.0;
        Vl << 0.0, 2.0;

        double vel[3] = { 1.0 + 1e-6, 2.0, 4.0 };
        double toi[3] = { 1.0, 0.5, 0.25 };
        double dx[3] = { 0.5, 0.0, -0.5 };

        SECTION("AllCases")
        {
            // touches, intersects, passes-trough
            for (int i = 0; i < 3; ++i) {
                // moving: both (same), ij, both (op), kl, both (same)
                // doesn't affect toi
                for (int j = 0; j < 5; ++j) {
                    // extension, no-deform, compression,
                    // doesn't affect toi
                    for (int k = 0; k < 3; ++k) {
                        Uk << 0.0, -(3 - j) * vel[i] / 2.0;
                        Ul << 0.0, -(3 - j) * vel[i] / 2.0;

                        Ui << -dx[k], (j - 1.0) * vel[i] / 2.0;
                        Uj << dx[k], (j - 1.0) * vel[i] / 2.0;

                        double toi_expected = toi[i];
                        std::string section = "vel=" + std::to_string(vel[i])
                            + " moving=" + std::to_string(j)
                            + " dx=" + std::to_string(dx[k]);

                        SECTION(section)
                        {
                            check_toi(Vi, Vj, Vk, Ui, Uj, Uk, toi_expected,
                                /*v2=*/true);
                        }

                        SECTION("flipped " + section)
                        {
                            // change order of edge indices (i.e edge symmetry)
                            check_toi(Vj, Vi, Vk, Uj, Ui, Uk, toi_expected,
                                /*v2=*/true);
                        }
                    }
                }
            }
        }

        SECTION("GradientVertexMoving")
        {
            vel[0] = 1.0
                + 1e-4; // otherwise finite diff gives no-impact on one side.

            for (int i = 0; i < 3; ++i) {
                Uk << 0.0, -vel[i];
                Ul << 0.0, -vel[i];

                Ui << 0.0, 0.0;
                Uj << 0.0, 0.0;

                SECTION("VertexMoving vel=" + std::to_string(vel[i]))
                {
                    check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
                }
            }
        }

        SECTION("GradientBothMoving")
        {
            // otherwise finite diff gives no-impact on one side.
            vel[0] = 1.0 + 1e-4;

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 5; ++j) {
                    Uk << 0.0, -(3 - j) * vel[i] / 2.0;
                    Ul << 0.0, -(3 - j) * vel[i] / 2.0;

                    Ui << 0.0, (j - 1) * vel[i] / 2.0;
                    Uj << 0.0, (j - 1) * vel[i] / 2.0;

                    SECTION("GradientBothMoving vel=" + std::to_string(vel[i]) + "j" + std::to_string(j))
                    {
                        check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
                    }
                }
            }
        }
        SECTION("GradientExtending"){
            Uk << 0.0, -2.0;
            Ul << 0.0, -2.0;

            for (int k = 0; k < 3; ++k) {
                Ui << -dx[k],0.0;
                Uj << dx[k], 0.0;

                SECTION("GradientExtending ex=" + std::to_string(dx[k]))
                {
                    check_gradient(Vi, Vj, Vk, Ui, Uj, Uk);
                }
            }
        }
    }

    SECTION("TangentImpact") //  (alpha=0 || alpha = 1)
    {

        Vi << -0.5, 0.0;
        Vj << -1.5, 0.0;
        Vk << 0.5, 0.0;
        Vl << 1.5, 0.0;

        double vel[3] = { 1.0, 2.0, 4.0 };
        double toi[3] = { 1.0, 0.5, 0.25 };
        double dx[3] = { 0.5, 0.0, -0.5 };

        // touches, intersects, passes-trough
        for (int i = 0; i < 3; ++i) {
            // moving: both (same), ij, both (op), kl, both (same)
            // doesn't affect toi
            for (int j = 0; j < 5; ++j) {
                for (int k = 0; k < 3; ++k) {
                    // extension, no-deform, compression,
                    // doesn't affect toi
                    SECTION("vel=" + std::to_string(vel[i]) + " moving="
                        + std::to_string(j) + " dx=" + std::to_string(dx[k]))
                    {
                        Uk << -(3 - j) * vel[i] / 2.0, 0.0;
                        Ul << -(3 - j) * vel[i] / 2.0, 0.0;

                        Ui << (j - 1.0) * vel[i] / 2.0, 0.0;
                        // we only move one so we don't change the toi
                        Uj << (j - 1.0) * vel[i] / 2.0, dx[k];

                        double toi_expected = toi[i];
                        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, toi_expected);

                        // change order of edge indices (i.e edge symmetry)
                        check_toi(Vj, Vi, Vk, Uj, Ui, Uk, toi_expected);
                    }
                }
            }
        }
    }
    SECTION("Double Impact") // (rotating edge)
    {
        // See fixtures/double-impact.json

        Vi << -1.0, 0.0;
        Vj << 1.0, 0.0;
        Vk << 0.0, 0.5;
        Vl << 0.0, 1.5;

        Ui << 1.6730970740318298, 0.8025388419628143;
        Uj << -1.616142749786377, -0.6420311331748962;
        Uk << 0.0, -1.0;
        Ul << 0.0, -1.0;

        check_toi(Vi, Vj, Vk, Ui, Uj, Uk, 0.4482900963, /*v2=*/true);
    }
}
