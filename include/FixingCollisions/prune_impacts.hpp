/**
    Detection collisions between different geometry.
    Includes continous collision detection to compute the time of impact.
    Supported geometry: point vs edge
*/

#ifndef PRUNE_IMPACTS_H
#define PRUNE_IMPACTS_H

#include <FixingCollisions/impact.hpp>

namespace ccd {

EdgeEdgeImpactsPtr prune_impacts(const EdgeEdgeImpactsPtr all_impacts);

}

#endif
