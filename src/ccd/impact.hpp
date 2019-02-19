/** Data structures for impacts between different geometry. */

#ifndef IMPACT_H
#define IMPACT_H

#include <Eigen/Core>
#include <vector>

namespace ccd {

// Forward declarations and typedefs
class Impact;
typedef std::shared_ptr<Impact> ImpactPtr;
typedef std::vector<ImpactPtr> Impacts;
typedef std::shared_ptr<Impacts> ImpactsPtr;

class EdgeVertexImpact;
typedef std::shared_ptr<EdgeVertexImpact> EdgeVertexImpactPtr;
typedef std::vector<EdgeVertexImpactPtr> EdgeVertexImpacts;
typedef std::shared_ptr<EdgeVertexImpacts> EdgeVertexImpactsPtr;

class EdgeEdgeImpact;
typedef std::shared_ptr<EdgeEdgeImpact> EdgeEdgeImpactPtr;
typedef std::vector<EdgeEdgeImpactPtr> EdgeEdgeImpacts;
typedef std::shared_ptr<EdgeEdgeImpacts> EdgeEdgeImpactsPtr;

/** Class representing an abstract impact. */
class Impact {
public:
    /** Time of impact. */
    const double time;

    /**
    Construct an abstract impact.

    @param time The time of impact.
    */
    Impact(double time);
    virtual ~Impact() = default;

    /**
    Compare two impacts to determine if impact0 comes before impact1.

    @param impact0 Impact to check if is came first.
    @param impact1 Impact to check if is came second.
    @return A boolean for if impact0->time <= impact1->time.
    */
    static bool compare_impacts_by_time(ImpactPtr impact0, ImpactPtr impact1);
};

/** Class representing an impact of an edge and vertex. */
class EdgeVertexImpact : public Impact {
public:
    /** Impacted edge. */
    const int edge_index;
    /** Parameter along the edge where the impact occured */
    const double alpha;
    /** Impacting vertex. */
    const int vertex_index;

    /**
    Construct an edge-vertex impact.

    @param time The time of impact.
    @param edge_index The index of the edge being impacted by a vertex.
    @param alpha The prameter along the edge where the impact occurred.
    @param vertex_index The index of the vertex impacting the edge.
    */
    EdgeVertexImpact(
        double time, int edge_index, double alpha, int vertex_index);

    /**
    Construct an edge-vertex impact from an edge-edge impact.

    @param edges A matrix where rows are edges and columns are vertex indices.
    @param ee_impact Edge-edge impact to convert to an edge-vertex impact.
    */
    EdgeVertexImpact(
        const Eigen::MatrixX2i& edges, const EdgeEdgeImpactPtr ee_impact);
};

/** Class representing an impact of an edge and edge. */
class EdgeEdgeImpact : public Impact {
public:
    /** Impacted edge. */
    const int impacted_edge_index;
    /** Parameter along the impacted edge where the impact occured. */
    const double impacted_alpha;
    /** Impacting edge. */
    const int impacting_edge_index;
    /** Parameter along the impacting edge where the impact occured. */
    const double impacting_alpha;

    /**
    Construct an edge-edge impact.

    @param time The time of impact.
    @param impacted_edge_index The index of the edge being impacted by another
        edge.
    @param impacted_alpha The parameter along the impacted edge where the impact
        occured.
    @param impacting_edge_index The index of the edge impacting the impacted
        edge.
    @param impacting_alpha The parameter along the impacting edge where the
        impact occured.
    */
    EdgeEdgeImpact(double time, int impacted_edge_index, double impacted_alpha,
        int impacting_edge_index, double impacting_alpha);

    /**
    Convert all edge-vertex impacts to correspoding edge-edge impacts. There may
    be multiple edge-edge impacts per edge-vertex impact depending on the
    connectivity.

    @param edges The matrix of edges where each row is two indices for the
        endpoints in the vertices matrix (not a prameter).
    @param ev_impacts Vector of edge-vertex impacts to convert to edge-edge
        impacts.
    @returns A EdgeEdgeImpactsPtr containing all of the edge-edge impacts.
    */
    static EdgeEdgeImpactsPtr convert_edge_vertex_to_edge_edge_impacts(
        const Eigen::MatrixX2i& edges, const EdgeVertexImpactsPtr ev_impacts);
};

}

#endif
