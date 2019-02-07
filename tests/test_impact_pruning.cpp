#include <iostream>
#include <stdlib.h> /* srand, rand */

#include <catch2/catch.hpp>

#include <FixingCollisions/impact.hpp>

using namespace ccd;

TEST_CASE("Test edge vertex impact to edge edge impact", "[impact_pruning]")
{
    // Create the impact manually
    double time = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    double alpha0 = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    int edge0_index = 0, edge1_index = 1, vertex_index = 2;
    EdgeVertexImpactsPtr ev_impacts(new EdgeVertexImpacts);
    ev_impacts->push_back(EdgeVertexImpactPtr(
        new EdgeVertexImpact(time, edge0_index, alpha0, vertex_index)));

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
