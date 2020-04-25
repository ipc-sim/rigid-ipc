// Test the mass utilities.

#include <catch2/catch.hpp>

#include <physics/mass.hpp>

TEST_CASE("Center of Mass", "[physics][mass]")
{
    using namespace ccd::physics;
    int dim = GENERATE(2);
    int num_vertices = 6 * GENERATE(10, 50, 100);

    // Test vertices positions for given rb position
    Eigen::MatrixXd vertices = Eigen::MatrixXd::Random(num_vertices, dim);

    Eigen::MatrixXi edges =
        Eigen::VectorXd::LinSpaced(num_vertices, 0, num_vertices - 1)
            .cast<int>();
    edges = Eigen::MatrixXi(
        Eigen::Map<Eigen::MatrixXi>(edges.data(), num_vertices / dim, dim));

    double total_mass1;
    Eigen::VectorX3d center_of_mass1;
    Eigen::MatrixXX3d moment_of_inertia1;
    compute_mass_properties(
        vertices, edges, total_mass1, center_of_mass1, moment_of_inertia1);

    double total_mass2 = compute_total_mass(vertices, edges);
    Eigen::VectorX3d center_of_mass2 = compute_center_of_mass(vertices, edges);
    Eigen::MatrixXX3d moment_of_inertia2 =
        compute_moment_of_inertia(vertices, edges);

    CAPTURE(dim, num_vertices);
    CHECK(total_mass1 == total_mass2);
    CHECK((center_of_mass1 - center_of_mass2).norm() < 1e-12);
    CHECK((moment_of_inertia1 - moment_of_inertia2).norm() < 1e-12);
}

// TODO: Add 3D RB test
