#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch.hpp>

#include <ccd/impact.hpp>
#include <ccd/prune_impacts.hpp>

using namespace ccd;

inline double rand_norm_double()
{
    return static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}

TEST_CASE("Test edge vertex impact to edge edge impact", "[impact_pruning]")
{
    // Create the impact manually
    double time = rand_norm_double(), impacted_alpha = rand_norm_double();
    int impacted_edge_index = 0, impacting_edge_index = 1, vertex_index = 2;
    EdgeVertexImpacts ev_impacts;
    ev_impacts.push_back(
        { time, impacted_edge_index, impacted_alpha, vertex_index });
    EdgeEdgeImpacts ee_impacts;
    SECTION("Free point contacting edge")
    {
        Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
        edges.row(impacted_edge_index) << 0, 1;
        convert_edge_vertex_to_edge_edge_impacts(edges, ev_impacts, ee_impacts);
        REQUIRE(ee_impacts.size() == 0);
    }
    SECTION("Only two edges contacting")
    {
        Eigen::Matrix<int, 2, 2, Eigen::RowMajor> edges;
        edges.row(impacted_edge_index) << 0, 1;

        SECTION("impacting_alpha = 0")
        {
            edges.row(impacting_edge_index) << vertex_index, 3;
            convert_edge_vertex_to_edge_edge_impacts(
                edges, ev_impacts, ee_impacts);
            REQUIRE(ee_impacts.size() == 1);
            CHECK(ee_impacts.front().impacting_alpha == Approx(0.0).margin(1e-3));
        }
        SECTION("impacting_alpha = 1")
        {
            edges.row(impacting_edge_index) << 3, vertex_index;
            convert_edge_vertex_to_edge_edge_impacts(
                edges, ev_impacts, ee_impacts);
            REQUIRE(ee_impacts.size() == 1);
            CHECK(ee_impacts.front().impacting_alpha == Approx(1.0));
        }
        CHECK(ee_impacts.front().time == Approx(time));
        CHECK(ee_impacts.front().impacted_edge_index == impacted_edge_index);
        CHECK(ee_impacts.front().impacted_alpha == Approx(impacted_alpha));
        CHECK(ee_impacts.front().impacting_edge_index == impacting_edge_index);
    }

    SECTION("Three edges contacting")
    {
        SECTION("Out of three edges")
        {
            Eigen::Matrix<int, 3, 2, Eigen::RowMajor> edges;
            edges.row(impacted_edge_index) << 0, 1;

            SECTION("impacting_alpha = 0")
            {
                edges.row(impacting_edge_index) << vertex_index, 3;
                edges.row(impacting_edge_index + 1) << vertex_index, 4;
                convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts, ee_impacts);
                CHECK(ee_impacts[0].impacting_alpha == Approx(0.0));
                CHECK(ee_impacts[1].impacting_alpha == Approx(0.0));
            }
            SECTION("impacting_alpha = 1")
            {
                edges.row(impacting_edge_index) << 3, vertex_index;
                edges.row(impacting_edge_index + 1) << 4, vertex_index;
                convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts, ee_impacts);
                CHECK(ee_impacts[0].impacting_alpha == Approx(1.0));
                CHECK(ee_impacts[1].impacting_alpha == Approx(1.0));
            }
        }
        SECTION("Out of four edges")
        {
            Eigen::Matrix<int, 4, 2, Eigen::RowMajor> edges;
            edges.row(impacted_edge_index) << 0, 1;
            edges.row(3) << 3, 4;

            SECTION("impacting_alpha = 0")
            {
                edges.row(impacting_edge_index) << vertex_index, 3;
                edges.row(impacting_edge_index + 1) << vertex_index, 4;
                convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts, ee_impacts);
                CHECK(ee_impacts[0].impacting_alpha == Approx(0));
                CHECK(ee_impacts[1].impacting_alpha == Approx(0));
            }
            SECTION("impacting_alpha = 1")
            {
                edges.row(impacting_edge_index) << 3, vertex_index;
                edges.row(impacting_edge_index + 1) << 4, vertex_index;
                convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts, ee_impacts);
                CHECK(ee_impacts[0].impacting_alpha == Approx(1));
                CHECK(ee_impacts[1].impacting_alpha == Approx(1));
            }
        }
        REQUIRE(ee_impacts.size() == 2);
        CHECK(ee_impacts[0].time == Approx(time));
        CHECK(ee_impacts[1].time == Approx(time));
        CHECK(ee_impacts[0].impacted_edge_index == impacted_edge_index);
        CHECK(ee_impacts[1].impacted_edge_index == impacted_edge_index);
        CHECK(ee_impacts[0].impacted_alpha == Approx(impacted_alpha));
        CHECK(ee_impacts[1].impacted_alpha == Approx(impacted_alpha));
        CHECK(ee_impacts[0].impacting_edge_index == impacting_edge_index);
        CHECK(ee_impacts[1].impacting_edge_index == impacting_edge_index + 1);
    }
}

TEST_CASE("Test impact pruning", "[impact_pruning]")
{
    // Create the impact manually
    double early_time = rand_norm_double() / 2,
           late_time = rand_norm_double() / 2 + 0.5;
    int impacted_edge_index = 0, impacting_edge_index = 1, edge2_index = 2;
    EdgeEdgeImpacts all_impacts;

    EdgeEdgeImpact early_impact = { early_time, impacted_edge_index,
        rand_norm_double(), impacting_edge_index, 0.0 };
    EdgeEdgeImpact late_impact = { late_time, impacted_edge_index,
        rand_norm_double(), edge2_index, 0.0 };
    EdgeToImpactMap pruned_impacts;
    SECTION("Two impacts against edge 0")
    {
        SECTION("Early then late impact")
        {
            all_impacts.push_back(early_impact);
            all_impacts.push_back(late_impact);
            prune_impacts(all_impacts, pruned_impacts);
        }
        SECTION("Late then early impact")
        {
            all_impacts.push_back(late_impact);
            all_impacts.push_back(early_impact);
            prune_impacts(all_impacts, pruned_impacts);
        }
        REQUIRE(pruned_impacts.size() == 3);
        for (auto impact : pruned_impacts) {
            auto other
                = impact.first == edge2_index ? late_impact : early_impact;
            CHECK(impact.second.time == Approx(other.time));
            CHECK(
                impact.second.impacted_edge_index == other.impacted_edge_index);
            CHECK(impact.second.impacted_alpha == Approx(other.impacted_alpha));
            CHECK(impact.second.impacting_edge_index
                == other.impacting_edge_index);
            CHECK(
                impact.second.impacting_alpha == Approx(other.impacting_alpha));
        }
    }
}
