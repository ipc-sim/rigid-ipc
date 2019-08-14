// Data structures for impacts between different geometry.
#include <ccd/impact.hpp>

namespace ccd {

EdgeVertexImpact::EdgeVertexImpact() {}

EdgeVertexImpact::EdgeVertexImpact(
    double time, long edge_index, double alpha, long vertex_index)
    : time(time)
    , edge_index(edge_index)
    , alpha(alpha)
    , vertex_index(vertex_index)
{
}

EdgeEdgeImpact::EdgeEdgeImpact() {}

EdgeEdgeImpact::EdgeEdgeImpact(double time, long impacted_edge_index,
    double impacted_alpha, long impacting_edge_index, double impacting_alpha)
    : time(time)
    , impacted_edge_index(impacted_edge_index)
    , impacted_alpha(impacted_alpha)
    , impacting_edge_index(impacting_edge_index)
    , impacting_alpha(impacting_alpha)
{
}

bool EdgeVertexImpact::operator==(const EdgeVertexImpact& other) const
{
    return this->time == other.time && this->edge_index == other.edge_index
        && this->alpha == other.alpha
        && this->vertex_index == other.vertex_index;
}

// Convert all edge-vertex impacts to correspoding edge-edge impacts. There may
// be multiple edge-edge impacts per edge-vertex impact depending on the
// connectivity.
void convert_edge_vertex_to_edge_edge_impacts(const Eigen::MatrixX2i& edges,
    const EdgeVertexImpacts& ev_impacts, EdgeEdgeImpacts& ee_impacts)
{
    ee_impacts.clear();
    for (EdgeVertexImpact ev_impact : ev_impacts) {
        for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
            if (edges(edge_index, 0) == ev_impact.vertex_index
                || edges(edge_index, 1) == ev_impact.vertex_index) {
                ee_impacts.push_back(EdgeEdgeImpact(ev_impact.time,
                    ev_impact.edge_index, ev_impact.alpha, edge_index,
                    edges(edge_index, 0) == ev_impact.vertex_index ? 0.0
                                                                   : 1.0));
            }
        }
    }
}

// Convert all edge-edge impacts to correspoding edge-vertex impacts. There may
// be multiple edge-vertex impacts per edge-edge impact depending on the
// connectivity.
void convert_edge_edge_to_edge_vertex_impacts(const Eigen::MatrixX2i& edges,
    const EdgeEdgeImpacts& ee_impacts, EdgeVertexImpacts& ev_impacts)
{
    ev_impacts.clear();
    for (EdgeEdgeImpact ee_impact : ee_impacts) {
        ev_impacts.push_back(EdgeVertexImpact(ee_impact.time,
            ee_impact.impacted_edge_index, ee_impact.impacted_alpha,
            edges(ee_impact.impacting_edge_index,
                int(ee_impact.impacting_alpha + 0.5))));
    }
}
} // namespace ccd
