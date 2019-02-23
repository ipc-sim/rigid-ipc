/// Prune impacts based on earliest contacts.
#ifndef PRUNE_IMPACTS_H
#define PRUNE_IMPACTS_H

#include <unordered_map>

#include <ccd/impact.hpp>

namespace ccd {

/// An unordered map from an edge index to an edge edge impact pointer.
typedef std::unordered_map<int, EdgeEdgeImpact> EdgeToImpactMap;

/// Prune the impacts to only include the earliest impact for each edge.
///
/// @param all_impacts All of the impacts to be pruned.
/// @param pruned_impacts Reference for the EdgeToImpactMap for where to store
///     the pruned impacts.
/// @return A subset of provided impacts where only the earliest impact per edge
///     is included are stored in pruned_impacts.
void prune_impacts(
    const EdgeEdgeImpacts& all_impacts, EdgeToImpactMap& pruned_impacts);

}

#endif
