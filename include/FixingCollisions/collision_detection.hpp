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

enum DetectionMethod { BRUTE_FORCE, HASH_MAP };

// Compute the time of impact of a point and edge moving in 2D.
double compute_edge_vertex_time_of_impact(const Eigen::MatrixX2d& vertex0,
    const Eigen::MatrixX2d& displacement0, const Eigen::MatrixX2d& edge_vertex1,
    const Eigen::MatrixX2d& edge_displacement1,
    const Eigen::MatrixX2d& edge_vertex2,
    const Eigen::MatrixX2d& edge_displacement2);

// General detection algorithm
ImpactsPtr detect_edge_vertex_collisions(const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXd& displacements, const Eigen::MatrixX2i& edges,
    DetectionMethod method = BRUTE_FORCE);

// Brute-force detection algorithm
ImpactsPtr detect_edge_vertex_collisions_brute_force(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges);

// Spatial-hashing detection algorithm
ImpactsPtr detect_edge_vertex_collisions_hash_map(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXd& displacements,
    const Eigen::MatrixX2i& edges);

} // namespace ccd
#endif
