/**
    Detection collisions between different geometry.
    Includes continous collision detection to compute the time of impact.
    Supported geometry: point vs edge
*/

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

    EdgeVertexImpact(
        double time, int edge_index, double alpha, int vertex_index);
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

    EdgeEdgeImpact(double time, int edge0_index, double alpha0, int edge1_index,
        double alpha1);

    static std::shared_ptr<std::vector<std::shared_ptr<EdgeEdgeImpact>>>
    convert_edge_vertex_to_edge_edge_impacts(const Eigen::MatrixX2i& edges,
        const EdgeVertexImpactsPtr other_impacts);
};

typedef std::shared_ptr<EdgeEdgeImpact> EdgeEdgeImpactPtr;
typedef std::vector<EdgeEdgeImpactPtr> EdgeEdgeImpacts;
typedef std::shared_ptr<EdgeEdgeImpacts> EdgeEdgeImpactsPtr;

}

#endif
