#include <iomanip>
#include <iostream>

#include <catch.hpp>

#include <autodiff/finitediff.hpp>
#include <autogen/collision_volume.hpp>
#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include "collision_generator.hpp"

using namespace ccd;

// ---------------------------------------------------
// Helper functions
// ---------------------------------------------------
/// Compares the volume of different implementations
/// against the expected volume.
/// Functions:
///     - ccd::collision_volume
///     - ccd::autogen::collision_volume
///     - ccd::autodiff::collision_volume<double>
///     - ccd::autodiff::collision_volume<DScalar>
/// The last two are higher-level functions and will
/// recompute toi & alpha.

void check_volume(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk, const double toi,
    const double alpha, const double epsilon, const double expected_vol);

/// Compares the autodiff gradient vs finite differences
bool check_gradient(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk, const double epsilon);

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("CollisionVolume", "[collision_volume][no-grad]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities
    double epsilon = 1.0;

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

        double edge_length = 2.0;
        double alpha = 0.5;

        SECTION("VertexMoving")
        {
            // touches, intersects, passes-trough
            auto vel = GENERATE(1.1, 2.0, 4.0);
            Uk << 0.0, -vel;
            double toi = 1.0 / vel;

            // volume = - (1.0 - t) * eps * || edge ||
            double expected_vol = -(1.0 - toi) * epsilon * edge_length;
            check_volume(
                Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha, epsilon, expected_vol);
        }

        SECTION("BothMoving")
        {
            // touches, intersects, passes-trough
            auto vel = GENERATE(1.1, 2.0, 4.0);
            // moving direction:
            //  -> ->, -> -, -> <-, - <-, <- <-
            auto j = GENERATE(0, 1, 2, 3, 4);

            Uk << 0.0, -(3 - j) * vel / 2.0;
            Ui << 0.0, (j - 1) * vel / 2.0;
            Uj << 0.0, (j - 1) * vel / 2.0;

            double toi = 1.0 / vel;

            // volume = - (1.0 - t) * sqrt( eps^2 + u_ij^2 ) * || edge ||
            double expected_vol = -(1.0 - toi)
                * std::sqrt(epsilon * epsilon + Ui(1) * Ui(1)) * edge_length;
            check_volume(
                Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha, epsilon, expected_vol);
        }
    }
}

TEST_CASE("CollisionVolumeGradient", "[collision_volume][gradient]")
{
    Eigen::Vector2d Vi, Vj, Vk; // positions
    Eigen::Vector2d Ui, Uj, Uk; // velocities
    double epsilon = 1.0;

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
            auto vel = GENERATE(1.1, 2.0, 4.0);
            Uk << 0.0, -vel;
            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
        }

        SECTION("BothMoving")
        {
            // touches, intersects, passes-trough
            auto vel = GENERATE(1.1, 2.0, 4.0);
            // moving direction:
            //  -> ->, -> -, -> <-, - <-, <- <-
            auto j = GENERATE(0, 1, 2, 3, 4);

            Uk << 0.0, -(3 - j) * vel / 2.0;

            Ui << 0.0, (j - 1) * vel / 2.0;
            Uj << 0.0, (j - 1) * vel / 2.0;

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
        }

        SECTION("Extending")
        {
            auto dx = GENERATE(0.5, 0.0, -0.5);
            Uk << 0.0, -2.0;

            Ui << -dx, 0.0;
            Uj << dx, 0.0;
            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
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

            check_gradient(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
        }
    }
}
TEST_CASE("CollisionVolumeRandom", "[collision_volume][random]")
{

    SECTION("GradientRigid_EPS1")
    {
        using namespace ccd::unittests;
        double epsilon = 1.0;
        auto impact = GENERATE(random_impacts(100, /*rigid=*/true));
        check_gradient(impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj,
            impact.Uk, epsilon);
    }

    SECTION("GradientNonRigid_EPS1")
    {
        using namespace ccd::unittests;
        double epsilon = 1.0;
        auto impact = GENERATE(random_impacts(100, /*rigid=*/false));

        if (!check_gradient(impact.Vi, impact.Vj, impact.Vk, impact.Ui, impact.Uj,
                            impact.Uk, epsilon)){
            std::cout << impact << std::endl;
        }
    }

}

void check_volume(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk, const double toi,
    const double alpha, const double epsilon, const double expected_vol)
{

    double actual_vol
        = ccd::collision_volume(Vi, Vj, Ui, Uj, toi, alpha, epsilon);
    CHECK(actual_vol == Approx(expected_vol));

    double autogen_vol
        = ccd::autogen::collision_volume(Vi, Vj, Ui, Uj, toi, alpha, epsilon);
    CHECK(autogen_vol == Approx(expected_vol));

    double autodiff_vol
        = ccd::autodiff::collision_volume(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
    CHECK(autodiff_vol == Approx(expected_vol));

    // check with autodiff variables
    DiffScalarBase::setVariableCount(8);

    DVector2 DUi = dvector(0, Ui);
    DVector2 DUj = dvector(2, Uj);
    DVector2 DUk = dvector(4, Uk);
    DScalar(6, 0.0);
    DScalar(7, 0.0);

    DScalar dautodiff_vol
        = ccd::autodiff::collision_volume(Vi, Vj, Vk, DUi, DUj, DUk, epsilon);
    CHECK(dautodiff_vol.getValue() == Approx(expected_vol));
}

bool check_gradient(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
    const Eigen::Vector2d& Vk, const Eigen::Vector2d& Ui,
    const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk, const double epsilon)
{
    Vector8d grad, grad_fd;
    grad
        = ccd::autodiff::collision_volume_grad(Vi, Vj, Vk, Ui, Uj, Uk, epsilon);
    REQUIRE(!std::isnan(grad.squaredNorm()));

    grad_fd = ccd::autodiff::collision_volume_grad_fd(
        Vi, Vj, Vk, Ui, Uj, Uk, epsilon);

    REQUIRE(!std::isnan(grad_fd.squaredNorm()));

    bool compare_grad = compare_gradient(grad, grad_fd);
    if (!compare_grad) {
        std::cout << " grad    " << grad.transpose() << std::endl;
        std::cout << " grad_fd " << grad_fd.transpose() << std::endl;
    }
    CHECK(compare_grad);
    return compare_grad;
}
