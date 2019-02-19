/** Prune impacts based on earliest contacts. */
#ifndef PRUNE_IMPACTS_H
#define PRUNE_IMPACTS_H

#include <unordered_map>

#include <ccd/impact.hpp>

namespace ccd {

/**
Prune the impacts to only include the earliest impact for each edge.

@param all_impacts All of the impacts to be pruned.
@return A subset of provided impacts where only the earliest impact per edge is
    included.
*/
std::shared_ptr<std::unordered_map<int, EdgeEdgeImpactPtr>> prune_impacts(
    const EdgeEdgeImpactsPtr all_impacts);

}

#endif
