#include <catch2/catch.hpp>

#include <state.hpp>

using namespace ccd;
using namespace opt;

TEST_CASE("hash grid on links", "[ccd][hashgrid]")
{
    State state;
    state.load_scene(std::string(FIXTURES_DIR) + "/chain/two-links.json");

    for (int i = 0; i < 10; i++) {
        state.displacements.setRandom();
        state.displacements *= 3;

        EdgeVertexImpacts brute_force_ev_impacts;
        detect_edge_vertex_collisions(state.vertices, state.displacements,
            state.edges, brute_force_ev_impacts, DetectionMethod::BRUTE_FORCE);
        EdgeVertexImpacts hash_ev_impacts;
        detect_edge_vertex_collisions(state.vertices, state.displacements,
            state.edges, hash_ev_impacts, DetectionMethod::HASH_GRID);

        REQUIRE(brute_force_ev_impacts.size() == hash_ev_impacts.size());
        for (const auto& impact : brute_force_ev_impacts) {
            CHECK(std::find(
                      hash_ev_impacts.begin(), hash_ev_impacts.end(), impact)
                != hash_ev_impacts.end());
        }
    }
}
