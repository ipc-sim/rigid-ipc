// Data structures for impacts between different geometry.
#include <ccd/impact.hpp>

namespace ccd {

bool EdgeVertexImpact::operator==(const EdgeVertexImpact& other) const
{
    return this->time == other.time && this->edge_index == other.edge_index
        && this->alpha == other.alpha
        && this->vertex_index == other.vertex_index;
}

// Compare two edge-vertex impacts to determine if impact0 comes before impact1.
template <typename Impact>
bool compare_impacts_by_time(Impact impact0, Impact impact1)
{
    return impact0.time <= impact1.time;
}

// Template instantiations
template bool compare_impacts_by_time<EdgeVertexImpact>(
    EdgeVertexImpact impact0, EdgeVertexImpact impact1);
template bool compare_impacts_by_time<EdgeEdgeImpact>(
    EdgeEdgeImpact impact0, EdgeEdgeImpact impact1);

// Convert all edge-vertex impacts to correspoding edge-edge impacts. There may
// be multiple edge-edge impacts per edge-vertex impact depending on the
// connectivity.
void convert_edge_vertex_to_edge_edge_impacts(const Eigen::MatrixX2i& edges,
    const EdgeVertexImpacts& ev_impacts, EdgeEdgeImpacts& ee_impacts)
{
    ee_impacts.clear();
    EdgeEdgeImpact ee_impact;
    for (EdgeVertexImpact ev_impact : ev_impacts) {
        for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
            auto edge = edges.row(edge_index);
            if (edge(0) == ev_impact.vertex_index
                || edge(1) == ev_impact.vertex_index) {

                ee_impact.time = ev_impact.time;
                ee_impact.impacted_edge_index = ev_impact.edge_index;
                ee_impact.impacted_alpha = ev_impact.alpha;
                ee_impact.impacting_edge_index = edge_index;
                ee_impact.impacting_alpha
                    = edge(0) == ev_impact.vertex_index ? 0.0 : 1.0;

                ee_impacts.push_back(ee_impact);
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
    EdgeVertexImpact ev_impact;
    for (EdgeEdgeImpact ee_impact : ee_impacts) {
        ev_impact.time = ee_impact.time;
        ev_impact.edge_index = ee_impact.impacted_edge_index;
        ev_impact.alpha = ee_impact.impacted_alpha;
        ev_impact.vertex_index = edges.row(ee_impact.impacting_edge_index)(
            int(ee_impact.impacting_alpha + 0.5));

        ev_impacts.push_back(ev_impact);
    }
}
} // namespace ccd
