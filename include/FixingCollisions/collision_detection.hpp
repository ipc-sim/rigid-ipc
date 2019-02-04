// Detection collisions between different geometry.
// Includes continous collision detection to compute the time of impact.
// Supported geometry: point vs edge
#ifndef COLLISION_DETECTION
#define COLLISION_DETECTION

#include <Eigen/Core>
#include <vector>

namespace ccd {
struct Impact {
    int vertex_index;
    int edge_index;
    double time;
};

typedef std::shared_ptr<Impact> ImpactPtr;
typedef std::vector<ImpactPtr> Impacts;
typedef std::shared_ptr<Impacts> ImpactsPtr;

enum DetectionMethod { BRUTE_FORCE,
    HASH_MAP };

ImpactPtr detect_edge_vertex_collision(
    const Eigen::MatrixXd& vertex0_t0,
    const Eigen::MatrixXd& vertex0_t1,
    const Eigen::MatrixXd& edge_vertex1_t0,
    const Eigen::MatrixXd& edge_vertex1_t1,
    const Eigen::MatrixXd& edge_vertex2_t0,
    const Eigen::MatrixXd& edge_vertex2_t1);

// General detection algorithm
ImpactsPtr detect_edge_vertex_collisions(
    const Eigen::MatrixXd& vertices_t0,
    const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixX2i& edges,
    DetectionMethod method = BRUTE_FORCE);

// Brute-force detection algorithm
ImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices_t0, const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixX2i& edges);

// Spatial-hashing detection algorithm
ImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& vertices_t0, const Eigen::MatrixXd& vertices_t1,
    const Eigen::MatrixX2i& edges);

} // namespace ccd
#endif
