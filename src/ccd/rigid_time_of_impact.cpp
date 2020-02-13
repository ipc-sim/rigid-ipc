#include "rigid_time_of_impact.hpp"

#include <ccd/interval_root_finder.hpp>
#include <geometry/distance.hpp>
#include <logger.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

/// Find time-of-impact between two rigid bodies
bool compute_edge_vertex_time_of_impact(
    const physics::RigidBody& bodyA,
    const physics::Pose<double>& displacementA, // Displacement of bodyA
    const size_t& vertex_id,                    // In bodyA
    const physics::RigidBody& bodyB,
    const physics::Pose<double>& displacementB, // Displacement of bodyB
    const size_t& edge_id,                      // In bodyB
    double& toi)
{
    // TODO: 3D
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);

    const physics::Pose<Interval> bodyA_pose_interval =
        bodyA.pose.cast<Interval>();
    const physics::Pose<Interval> bodyB_pose_interval =
        bodyB.pose.cast<Interval>();

    const physics::Pose<Interval> displacementA_interval =
        displacementA.cast<Interval>();
    const physics::Pose<Interval> displacementB_interval =
        displacementB.cast<Interval>();

    const auto distance = [&](Interval t) {
        // Compute the poses at time t
        physics::Pose<Interval> poseA =
            bodyA_pose_interval + displacementA_interval * t;
        physics::Pose<Interval> poseB =
            bodyB_pose_interval + displacementB_interval * t;

        // Get the world vertex of the point at time t
        Eigen::VectorX3<Interval> v0 =
            bodyA.world_vertex<Interval>(poseA, vertex_id);
        // spdlog::debug("{}", logger::fmt_eigen_intervals(v0));
        // Get the world vertex of the edge at time t
        Eigen::VectorX3<Interval> v1 =
            bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 0));
        Eigen::VectorX3<Interval> v2 =
            bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 1));

        return geometry::point_line_signed_distance<Interval>(v0, v1, v2);
    };

    const auto is_point_along_edge = [&](Interval t) {
        // Compute the poses at time t
        physics::Pose<Interval> poseA =
            bodyA_pose_interval + displacementA_interval * t;
        physics::Pose<Interval> poseB =
            bodyB_pose_interval + displacementB_interval * t;

        // Get the world vertex of the point at time t
        Eigen::VectorX3<Interval> v0 =
            bodyA.world_vertex<Interval>(poseA, vertex_id);
        // Get the world vertex of the edge at time t
        Eigen::VectorX3<Interval> v1 =
            bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 0));
        Eigen::VectorX3<Interval> v2 =
            bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 1));

        // Project the point onto the edge by computing its scalar projection
        Eigen::VectorX3<Interval> edge_vec = v2 - v1;
        Interval alpha = (v0 - v1).dot(edge_vec) / edge_vec.squaredNorm();
        // spdlog::debug("α ∈ {}", logger::fmt_interval(alpha));
        return overlap(alpha, Interval(0, 1));
    };

    Interval toi_interval;
    // TODO: Set tolerance dynamically
    bool is_impacting = interval_root_finder(
        distance, is_point_along_edge, Interval(0, 1), toi_interval);
    toi = median(toi_interval);
    return is_impacting;
}

} // namespace ccd
