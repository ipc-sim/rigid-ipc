/// Prune impacts based on earliest contacts.
#ifndef PRUNE_IMPACTS_H
#define PRUNE_IMPACTS_H

#include <unordered_map>

#include <ccd/impact.hpp>

namespace ccd {

/// Prune the impacts to only include the earliest impact for each edge.
///
/// @param all_impacts All of the impacts to be pruned.
/// @param pruned_impact_indices A vector of lengeth the numbe of edges. Store
/// the index of the earliest impact in this vector.
/// @return A subset of provided impacts where only the earliest impact per edge
///     is included are stored in pruned_impact_indices.
void prune_impacts(
    const EdgeEdgeImpacts& all_impacts, Eigen::VectorXi& pruned_impact_indices);

}

#endif
