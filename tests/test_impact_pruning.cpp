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
    double time = rand_norm_double(), alpha0 = rand_norm_double();
    int edge0_index = 0, edge1_index = 1, vertex_index = 2;
    EdgeVertexImpactsPtr ev_impacts(new EdgeVertexImpacts);
    ev_impacts->push_back(EdgeVertexImpactPtr(
        new EdgeVertexImpact(time, edge0_index, alpha0, vertex_index)));

    SECTION("Free point contacting edge")
    {
        Eigen::Matrix<int, 1, 2, Eigen::RowMajor> edges;
        edges.row(edge0_index) << 0, 1;
        EdgeEdgeImpactsPtr ee_impacts
            = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                edges, ev_impacts);
        REQUIRE(ee_impacts->size() == 0);
    }
    SECTION("Only two edges contacting")
    {
        Eigen::Matrix<int, 2, 2, Eigen::RowMajor> edges;
        edges.row(edge0_index) << 0, 1;

        EdgeEdgeImpactsPtr ee_impacts;
        SECTION("alpha1 = 0")
        {
            edges.row(edge1_index) << vertex_index, 3;
            ee_impacts
                = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts);
            REQUIRE(ee_impacts->front()->alpha1 == Approx(0));
        }
        SECTION("alpha1 = 1")
        {
            edges.row(edge1_index) << 3, vertex_index;
            ee_impacts
                = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                    edges, ev_impacts);
            REQUIRE(ee_impacts->front()->alpha1 == Approx(1));
        }
        REQUIRE(ee_impacts->size() == 1);
        REQUIRE(ee_impacts->front()->time == Approx(time));
        REQUIRE(ee_impacts->front()->edge0_index == edge0_index);
        REQUIRE(ee_impacts->front()->alpha0 == Approx(alpha0));
        REQUIRE(ee_impacts->front()->edge1_index == edge1_index);
    }

    SECTION("Three edges contacting")
    {
        EdgeEdgeImpactsPtr ee_impacts;
        SECTION("Out of three edges")
        {
            Eigen::Matrix<int, 3, 2, Eigen::RowMajor> edges;
            edges.row(edge0_index) << 0, 1;

            SECTION("alpha1 = 0")
            {
                edges.row(edge1_index) << vertex_index, 3;
                edges.row(edge1_index + 1) << vertex_index, 4;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->alpha1 == Approx(0));
                REQUIRE((*ee_impacts)[1]->alpha1 == Approx(0));
            }
            SECTION("alpha1 = 1")
            {
                edges.row(edge1_index) << 3, vertex_index;
                edges.row(edge1_index + 1) << 4, vertex_index;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->alpha1 == Approx(1));
                REQUIRE((*ee_impacts)[1]->alpha1 == Approx(1));
            }
        }
        SECTION("Out of four edges")
        {
            Eigen::Matrix<int, 4, 2, Eigen::RowMajor> edges;
            edges.row(edge0_index) << 0, 1;
            edges.row(3) << 3, 4;

            SECTION("alpha1 = 0")
            {
                edges.row(edge1_index) << vertex_index, 3;
                edges.row(edge1_index + 1) << vertex_index, 4;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->alpha1 == Approx(0));
                REQUIRE((*ee_impacts)[1]->alpha1 == Approx(0));
            }
            SECTION("alpha1 = 1")
            {
                edges.row(edge1_index) << 3, vertex_index;
                edges.row(edge1_index + 1) << 4, vertex_index;
                ee_impacts
                    = EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
                        edges, ev_impacts);
                REQUIRE((*ee_impacts)[0]->alpha1 == Approx(1));
                REQUIRE((*ee_impacts)[1]->alpha1 == Approx(1));
            }
        }
        REQUIRE(ee_impacts->size() == 2);
        REQUIRE((*ee_impacts)[0]->time == Approx(time));
        REQUIRE((*ee_impacts)[1]->time == Approx(time));
        REQUIRE((*ee_impacts)[0]->edge0_index == edge0_index);
        REQUIRE((*ee_impacts)[1]->edge0_index == edge0_index);
        REQUIRE((*ee_impacts)[0]->alpha0 == Approx(alpha0));
        REQUIRE((*ee_impacts)[1]->alpha0 == Approx(alpha0));
        REQUIRE((*ee_impacts)[0]->edge1_index == edge1_index);
        REQUIRE((*ee_impacts)[1]->edge1_index == edge1_index + 1);
    }
}

TEST_CASE("Test impact pruning", "[impact_pruning]")
{
    // Create the impact manually
    double early_time = rand_norm_double() / 2,
           late_time = rand_norm_double() / 2 + 0.5;
    int edge0_index = 0, edge1_index = 1, edge2_index = 2;
    EdgeEdgeImpactsPtr all_impacts(new EdgeEdgeImpacts);

    EdgeEdgeImpactPtr early_impact(new EdgeEdgeImpact(
        early_time, edge0_index, rand_norm_double(), edge1_index, 0)),
        late_impact(new EdgeEdgeImpact(
            late_time, edge0_index, rand_norm_double(), edge2_index, 0));
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
