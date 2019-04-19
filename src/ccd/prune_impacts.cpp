// Prune impacts based on earliest contacts.
#include <ccd/prune_impacts.hpp>

namespace ccd {

// Prune the impacts so only the earliest impacts per edge are included.
int prune_impacts(
    const EdgeEdgeImpacts& all_impacts, Eigen::VectorXi& pruned_impact_indices)
{
    // An index value of -1 indicates no impact for that edge
    pruned_impact_indices.setConstant(-1);

    // Loop over all impacts and add them to the prunced impacts if they are the
    // earliest.
    EdgeEdgeImpact ee_impact;
    int num_pruned_impacts = 0;
    for (size_t i = 0; i < all_impacts.size(); i++) {
        ee_impact = all_impacts[i];
        for (int index :
            { ee_impact.impacted_edge_index, ee_impact.impacting_edge_index }) {
            // If the edge is not in the map or the current impact happens
            // later replace it in the pruned impacts.
            if (pruned_impact_indices[index] == -1
                || all_impacts[size_t(pruned_impact_indices[index])].time
                    > ee_impact.time) {
                num_pruned_impacts += int(pruned_impact_indices[index] == -1);
                pruned_impact_indices[index] = int(i);
            }
        }
    }

    return num_pruned_impacts;
}

} // namespace ccd
