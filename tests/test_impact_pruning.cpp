#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch.hpp>

#include <FixingCollisions/impact.hpp>
#include <FixingCollisions/prune_impacts.hpp>

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
    EdgeVertexImpactsPtr ev_impacts = std::make_shared<EdgeVertexImpacts>();
    ev_impacts->push_back(std::make_shared<EdgeVertexImpact>(
        time, impacted_edge_index, impacted_alpha, vertex_index));

    SECTION("Free point contacting edge")
    {
        Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
        edges.row(impacted_edge_index) << 0, 1;
        EdgeEdgeImpactsPtr ee_impacts
            = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                edges, ev_impacts);
        REQUIRE(ee_impacts->size() == 0);
    }
    SECTION("Only two edges contacting")
    {
        Eigen::Matrix<int, 2, 2, Eigen::RowMajor> edges;
        edges.row(impacted_edge_index) << 0, 1;

        EdgeEdgeImpactsPtr ee_impacts;
        SECTION("impacting_alpha = 0")
        {
            edges.row(impacting_edge_index) << vertex_index, 3;
            ee_impacts
                = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts);
            REQUIRE(ee_impacts->front()->impacting_alpha == Approx(0));
        }
        SECTION("impacting_alpha = 1")
        {
            edges.row(impacting_edge_index) << 3, vertex_index;
            ee_impacts
                = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts);
            REQUIRE(ee_impacts->front()->impacting_alpha == Approx(1));
        }
        REQUIRE(ee_impacts->size() == 1);
        REQUIRE(ee_impacts->front()->time == Approx(time));
        REQUIRE(
            ee_impacts->front()->impacted_edge_index == impacted_edge_index);
        REQUIRE(ee_impacts->front()->impacted_alpha == Approx(impacted_alpha));
        REQUIRE(
            ee_impacts->front()->impacting_edge_index == impacting_edge_index);
    }

    SECTION("Three edges contacting")
    {
        EdgeEdgeImpactsPtr ee_impacts;
        SECTION("Out of three edges")
        {
            Eigen::Matrix<int, 3, 2, Eigen::RowMajor> edges;
            edges.row(impacted_edge_index) << 0, 1;

            SECTION("impacting_alpha = 0")
            {
                edges.row(impacting_edge_index) << vertex_index, 3;
                edges.row(impacting_edge_index + 1) << vertex_index, 4;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->impacting_alpha == Approx(0));
                REQUIRE((*ee_impacts)[1]->impacting_alpha == Approx(0));
            }
            SECTION("impacting_alpha = 1")
            {
                edges.row(impacting_edge_index) << 3, vertex_index;
                edges.row(impacting_edge_index + 1) << 4, vertex_index;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->impacting_alpha == Approx(1));
                REQUIRE((*ee_impacts)[1]->impacting_alpha == Approx(1));
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
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->impacting_alpha == Approx(0));
                REQUIRE((*ee_impacts)[1]->impacting_alpha == Approx(0));
            }
            SECTION("impacting_alpha = 1")
            {
                edges.row(impacting_edge_index) << 3, vertex_index;
                edges.row(impacting_edge_index + 1) << 4, vertex_index;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->impacting_alpha == Approx(1));
                REQUIRE((*ee_impacts)[1]->impacting_alpha == Approx(1));
            }
        }
        REQUIRE(ee_impacts->size() == 2);
        REQUIRE((*ee_impacts)[0]->time == Approx(time));
        REQUIRE((*ee_impacts)[1]->time == Approx(time));
        REQUIRE((*ee_impacts)[0]->impacted_edge_index == impacted_edge_index);
        REQUIRE((*ee_impacts)[1]->impacted_edge_index == impacted_edge_index);
        REQUIRE((*ee_impacts)[0]->impacted_alpha == Approx(impacted_alpha));
        REQUIRE((*ee_impacts)[1]->impacted_alpha == Approx(impacted_alpha));
        REQUIRE((*ee_impacts)[0]->impacting_edge_index == impacting_edge_index);
        REQUIRE(
            (*ee_impacts)[1]->impacting_edge_index == impacting_edge_index + 1);
    }
}

TEST_CASE("Test impact pruning", "[impact_pruning]")
{
    // Create the impact manually
    double early_time = rand_norm_double() / 2,
           late_time = rand_norm_double() / 2 + 0.5;
    int impacted_edge_index = 0, impacting_edge_index = 1, edge2_index = 2;
    EdgeEdgeImpactsPtr all_impacts = std::make_shared<EdgeEdgeImpacts>();

    EdgeEdgeImpactPtr early_impact
        = std::make_shared<EdgeEdgeImpact>(early_time, impacted_edge_index,
            rand_norm_double(), impacting_edge_index, 0);
    EdgeEdgeImpactPtr late_impact = std::make_shared<EdgeEdgeImpact>(
        late_time, impacted_edge_index, rand_norm_double(), edge2_index, 0);
    SECTION("Two impacts against edge 0")
    {
        std::shared_ptr<std::unordered_map<int, EdgeEdgeImpactPtr>>
            pruned_impacts;
        SECTION("Early then late impact")
        {
            all_impacts->push_back(early_impact);
            all_impacts->push_back(late_impact);
            pruned_impacts = prune_impacts(all_impacts);
        }
        SECTION("Late then early impact")
        {
            all_impacts->push_back(late_impact);
            all_impacts->push_back(early_impact);
            pruned_impacts = prune_impacts(all_impacts);
        }
        REQUIRE(pruned_impacts->size() == 3);
        for (auto impact : *pruned_impacts) {
            if (impact.first == edge2_index) {
                REQUIRE(impact.second == late_impact);
            } else {
                REQUIRE(impact.second == early_impact);
            }
        }
    }
}
