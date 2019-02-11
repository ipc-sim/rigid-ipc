/** Data structures for impacts between different geometry. */

#ifndef IMPACT_H
#define IMPACT_H

#include <Eigen/Core>
#include <vector>

namespace ccd {

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
    virtual ~Impact();
};

typedef std::shared_ptr<Impact> ImpactPtr;
typedef std::vector<ImpactPtr> Impacts;
typedef std::shared_ptr<Impacts> ImpactsPtr;

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
    Compare two edge-vertex impacts to determine if impact0 comes before
    impact1.

    @param impact0 Impact to check if is came first.
    @param impact1 Impact to check if is came second.
    @return A boolean for if impact0->time <= impact1->time.
    */
    static bool compare_impacts_by_time(
        std::shared_ptr<EdgeVertexImpact> impact0,
        std::shared_ptr<EdgeVertexImpact> impact1);
};

typedef std::shared_ptr<EdgeVertexImpact> EdgeVertexImpactPtr;
typedef std::vector<EdgeVertexImpactPtr> EdgeVertexImpacts;
typedef std::shared_ptr<EdgeVertexImpacts> EdgeVertexImpactsPtr;

/** Class representing an impact of an edge and edge. */
class EdgeEdgeImpact : public Impact {
public:
    /** Impacting edge 0. */
    const int edge0_index;
    /** Parameter along edge0 where the impact occured */
    const double alpha0;
    /** Impacted edge 1. */
    const int edge1_index;
    /** Parameter along edge1 where the impact occured */
    const double alpha1;

    /**
    Construct an edge-edge impact.

    @param time The time of impact.
    @param edge0_index The index of the edge being impacted by another edge.
    @param alpha0 The prameter along edge0 where the impact occurred.
    @param edge1_index The index of the edge impacting the edge0.
    @param alpha1 The prameter along edge1 where the impact occurred.
    */
    EdgeEdgeImpact(double time, int edge0_index, double alpha0, int edge1_index,
        double alpha1);

    /**
    Convert all edge-vertex impacts to correspoding edge-edge impacts.
    There may be multiple edge-edge impacts per edge-vertex impact
    depending on the connectivity.

    @param edges The matrix of edges where each row is two indices for the
                 endpoints in the vertices matrix (not a prameter).
    @param ev_impacts Vector of edge-vertex impacts to convert to edge-edge
                      impacts.
    @returns A EdgeEdgeImpactsPtr containing all of the edge-edge impacts.
    */
    static std::shared_ptr<std::vector<std::shared_ptr<EdgeEdgeImpact>>>
    convert_edge_vertex_to_edge_edge_impacts(const Eigen::MatrixX2i& edges,
        const EdgeVertexImpactsPtr other_impacts);
};

typedef std::shared_ptr<EdgeEdgeImpact> EdgeEdgeImpactPtr;
typedef std::vector<EdgeEdgeImpactPtr> EdgeEdgeImpacts;
typedef std::shared_ptr<EdgeEdgeImpacts> EdgeEdgeImpactsPtr;

}

#endif
