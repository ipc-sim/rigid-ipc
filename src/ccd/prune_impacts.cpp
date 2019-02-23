// Prune impacts based on earliest contacts.
#include <ccd/prune_impacts.hpp>

namespace ccd {

// Prune the impacts so only the earliest impacts per edge are included.
void prune_impacts(
    const EdgeEdgeImpacts& all_impacts, EdgeToImpactMap& pruned_impacts)
{
    // Store ealiest impacts as a hashmap of edge index to impact.
    pruned_impacts.clear();

    // Loop over all impacts and add them to the prunced impacts if they are the
    // earliest.
    for (EdgeEdgeImpact ee_impact : all_impacts) {
        for (int index :
            { ee_impact.impacted_edge_index, ee_impact.impacting_edge_index }) {
            // Find the impact for the edge.
            auto search = pruned_impacts.find(index);
            // If the edge is not in the map or the current impact happens
            // later replace it in the pruned impacts.
            if (search == pruned_impacts.end()
                || search->second.time > ee_impact.time) {
                pruned_impacts[index] = ee_impact;
            }
        }
    }
}

}
