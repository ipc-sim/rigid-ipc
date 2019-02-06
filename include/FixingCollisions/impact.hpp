/**
    Detection collisions between different geometry.
    Includes continous collision detection to compute the time of impact.
    Supported geometry: point vs edge
*/

#ifndef IMPACT_H
#define IMPACT_H

#include <vector>

namespace ccd {

/** Structure representing an impact of an edge and vertex. */
class Impact {
public:
    /** Impacting vertex. */
    const int vertex_index;
    /** Impacted edge. */
    const int edge_index;
    /** Time of impact. */
    const double time;
    /** Parameter along the edge where the impact occured */
    const double alpha;

    Impact(int vertex_index, int edge_index, double time, double alpha);
    virtual ~Impact();
};

typedef std::shared_ptr<Impact> ImpactPtr;
typedef std::vector<ImpactPtr> Impacts;
typedef std::shared_ptr<Impacts> ImpactsPtr;

}

#endif
