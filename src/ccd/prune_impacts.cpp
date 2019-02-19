// Prune impacts based on earliest contacts.
#include <ccd/prune_impacts.hpp>

namespace ccd {

// Prune the impacts so only the earliest impacts per edge are included.
std::shared_ptr<std::unordered_map<int, EdgeEdgeImpactPtr>> prune_impacts(
    const EdgeEdgeImpactsPtr all_impacts)
{
    // Store ealiest impacts as a hashmap of edge index to impact.
    std::shared_ptr<std::unordered_map<int, EdgeEdgeImpactPtr>> pruned_impacts
        = std::make_shared<std::unordered_map<int, EdgeEdgeImpactPtr>>();

    // Loop over all impacts and add them to the prunced impacts if they are the
    // earliest.
    for (EdgeEdgeImpactPtr& impact : *all_impacts) {
        for (int index :
            { impact->impacted_edge_index, impact->impacting_edge_index }) {
            // Find the impact for the edge.
            auto search = pruned_impacts->find(index);
            // If the edge is not in the map or the current impact happens
            // later replace it in the pruned impacts.
            if (search == pruned_impacts->end()
                || search->second->time > impact->time) {
                (*pruned_impacts)[index] = impact;
            }
        }
    }

    return pruned_impacts;
    // Convert the pruned impacts to a EdgeEdgeImpactsPtr
    // EdgeEdgeImpactsPtr resulting_impacts =
    // std::make_shared<EdgeEdgeImpacts>();
    // for (auto pruned_impact : pruned_impacts) {
    //     resulting_impacts->push_back(pruned_impact.second);
    // }
    // return resulting_impacts;
}

}
