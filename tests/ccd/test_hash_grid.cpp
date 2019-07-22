#include <catch2/catch.hpp>

#include <fstream>
#include <iostream>

#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>
#include <nlohmann/json.hpp>

using namespace ccd;
using namespace nlohmann;

TEST_CASE("hash grid gets same impacts as brute force", "ccdhashgrid")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd displacements;
    vertices.resize(10, 2);
    vertices << -2.4375, 0.0, 0.0, 0.0, 2.4375, 0.0, -3.0, -3.0, 0.0, -3.0, 3.0,
        -3.0, -3.0, -6.0, -0.5625, -6.0, 0.5625, -6.0, 3.0, -6.0;

    edges.resize(9, 2);
    edges << 0, 1, 1, 2, 1, 4, 3, 4, 4, 5, 3, 6, 5, 9, 6, 7, 8, 9;

    displacements.resize(vertices.rows(), vertices.cols());
    for (int i = 0; i < 10; i++) {
        displacements.setRandom();
        displacements *= 3;

        EdgeVertexImpacts brute_force_ev_impacts;
        detect_edge_vertex_collisions(vertices, displacements, edges,
            brute_force_ev_impacts, DetectionMethod::BRUTE_FORCE);
        EdgeVertexImpacts hash_ev_impacts;
        detect_edge_vertex_collisions(vertices, displacements, edges,
            hash_ev_impacts, DetectionMethod::HASH_GRID);

        REQUIRE(brute_force_ev_impacts.size() == hash_ev_impacts.size());
        for (const auto& impact : brute_force_ev_impacts) {
            CHECK(std::find(
                      hash_ev_impacts.begin(), hash_ev_impacts.end(), impact)
                != hash_ev_impacts.end());
        }
    }
}
