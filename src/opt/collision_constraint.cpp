#include "collision_constraint.hpp"

#include <ccd/prune_impacts.hpp>

namespace ccd {
namespace opt {

    CollisionConstraint::CollisionConstraint()
        : detection_method(BRUTE_FORCE)
        , recompute_collision_set(true)
    {
    }

    CollisionConstraint::~CollisionConstraint() {}


    void CollisionConstraint::detecteCollisions(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk)
    {
        ccd::detect_edge_vertex_collisions(
            vertices, Uk, edges, ev_impacts, detection_method, false);
        ccd::convert_edge_vertex_to_edge_edge_impacts(edges, ev_impacts, ee_impacts);
        num_pruned_impacts = prune_impacts(ee_impacts, edge_impact_map);
    }

} // namespace opt
} // namespace ccd
