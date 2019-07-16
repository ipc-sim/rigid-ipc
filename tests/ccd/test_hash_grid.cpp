#include <catch2/catch.hpp>

#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>
#include <ccd/collision_detection.hpp>
#include <io/serialize_json.hpp>

using namespace ccd;
using namespace nlohmann;

TEST_CASE("hash grid gets same impacts as brute force", "[ccd][hashgrid]")
{
    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd displacements;

    std::ifstream input(
        std::string(FIXTURES_DIR) + "/chain/two-links-mesh.json");
    json scene = json::parse(input);

    ccd::io::from_json(scene["vertices"], vertices);
    ccd::io::from_json(scene["edges"], edges);

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
