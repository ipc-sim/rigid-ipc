// Prune impacts based on earliest contacts.
#include <ccd/prune_impacts.hpp>

namespace ccd {

// Prune the impacts so only the earliest impacts per edge are included.
void prune_impacts(
    const EdgeEdgeImpacts& all_impacts, Eigen::VectorXi& pruned_impact_indices)
{
    // An index value of -1 indicates no impact for that edge
    pruned_impact_indices.setConstant(-1);

    // Loop over all impacts and add them to the prunced impacts if they are the
    // earliest.
    EdgeEdgeImpact ee_impact;
    for (int i = 0; i < all_impacts.size(); i++) {
        ee_impact = all_impacts[i];
        for (int index :
            { ee_impact.impacted_edge_index, ee_impact.impacting_edge_index }) {
            // If the edge is not in the map or the current impact happens
            // later replace it in the pruned impacts.
            if (pruned_impact_indices[index] == -1
                || all_impacts[pruned_impact_indices[index]].time
                    > ee_impact.time) {
                pruned_impact_indices[index] = i;
            }
        }
    }
}

}
