#include <FixingCollisions/impact.hpp>

namespace ccd {

Impact::Impact(double time)
    : time(time)
{
}

Impact::~Impact() {}

EdgeVertexImpact::EdgeVertexImpact(
    double time, int edge_index, double alpha, int vertex_index)
    : Impact(time)
    , edge_index(edge_index)
    , alpha(alpha)
    , vertex_index(vertex_index)
{
}

EdgeEdgeImpact::EdgeEdgeImpact(
    double time, int edge0_index, double alpha0, int edge1_index, double alpha1)
    : Impact(time)
    , edge0_index(edge0_index)
    , alpha0(alpha0)
    , edge1_index(edge1_index)
    , alpha1(alpha1)
{
}

EdgeEdgeImpactsPtr EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
    const Eigen::MatrixX2i& edges, const EdgeVertexImpactsPtr ev_impacts)
{
    EdgeEdgeImpactsPtr ee_impacts(new EdgeEdgeImpacts());
    for (auto ev_impact : *ev_impacts) {
        for (int edge_index = 0; edge_index < edges.rows(); edge_index++) {
            auto edge = edges.row(edge_index);
            if (edge(0) == ev_impact->vertex_index
                || edge(1) == ev_impact->vertex_index) {
                ee_impacts->push_back(EdgeEdgeImpactPtr(new EdgeEdgeImpact(
                    ev_impact->time, ev_impact->edge_index, ev_impact->alpha,
                    edge_index, edge(0) == ev_impact->vertex_index ? 0 : 1)));
            }
        }
    }
    return ee_impacts;
}

}
