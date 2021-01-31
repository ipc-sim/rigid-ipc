#include <catch2/catch.hpp>

#include <ccd/ccd.hpp>
#include <ccd/linear/edge_vertex_ccd.hpp>

using namespace ipc::rigid;

TEST_CASE("Test Continous Collision Detection", "[ccd]")
{
    Eigen::Vector2d v_t0, v_t1, e0_t0, e0_t1, e1_t0, e1_t1;
    double toi_expected, alpha_expected;
    bool is_collision_expected;
    SECTION("Edge becomes degenerate")
    {
        v_t0 << 0, 1;
        e0_t0 << -1, 0;
        e1_t0 << 1, 0;
        SECTION("Edge degenerates before impact")
        {
            v_t1 << 0, -1;
            e0_t1 << 3, 0;
            e1_t1 << -3, 0;
            // The edge will become degenerate at t=0.25
        }
        SECTION("Edge degenerates after impact")
        {
            v_t1 << 0, -1;
            e0_t1 << 0.5, 0;
            e1_t1 << -0.5, 0;
            // The edge will become degenerate at t=2/3
        }
        // The point will collide with the edge at t=0.5
        is_collision_expected = true;
        toi_expected = 0.5;
        alpha_expected = 0.5;
    }
    SECTION("Edge moving right; point moving left")
    {
        v_t0 << -1, 0;
        e0_t0 << 1, -1;
        e1_t0 << 1, 1;

        v_t1 << 1, 0;
        e0_t1 << -1, -1;
        e1_t1 << -1, 1;

        is_collision_expected = true;
        toi_expected = 0.5;
        alpha_expected = 0.5;
    }
    SECTION("Point on edge's line moving towards edge")
    {
        v_t0 << 0, 0;
        e0_t0 << 0, 1;
        e1_t0 << 0, 2;

        v_t1 << 0, 2;
        e0_t1 << 0, 1;
        e1_t1 << 0, 2;

        is_collision_expected = true;
        toi_expected = 0.5;
        alpha_expected = 0;
    }
    SECTION("Point and edge moving parallel")
    {
        v_t0 << 0, 1;
        e0_t0 << 1, 0;
        e1_t0 << 1, 2;

        v_t1 << 0, 2;
        e0_t1 << 1, 1;
        e1_t1 << 1, 3;

        is_collision_expected = false;
    }
    SECTION("Point moving right; edge stretching vertically")
    {
        e1_t0 << 1, -1;
        e1_t1 << 1, -2;
        SECTION("Swap vertices order e_0 = [0, 2]")
        {
            v_t0 << 1, 1;
            e0_t0 << 0, 0;

            v_t1 << 1, 2;
            e0_t1 << 1, 0;

            std::swap(v_t0, e0_t0);
            std::swap(v_t1, e0_t1);
        }
        SECTION("Swap vertices order e_0 = [1, 2]")
        {
            v_t0 << 0, 0;
            e0_t0 << 1, 1;

            v_t1 << 1, 0;
            e0_t1 << 1, 2;
        }
        is_collision_expected = true;
        toi_expected = 1.0;
        alpha_expected = 0.5;
    }

    double toi, alpha;
    bool is_colliding = compute_edge_vertex_time_of_impact(
        e0_t0, e1_t0, v_t0, (e0_t1 - e0_t0).eval(), (e1_t1 - e1_t0).eval(),
        (v_t1 - v_t0).eval(), toi, alpha);
    REQUIRE(is_colliding == is_collision_expected);
    if (is_collision_expected) {
        CHECK(toi == Approx(toi_expected));
        CHECK(alpha == Approx(alpha_expected));
    }
}
