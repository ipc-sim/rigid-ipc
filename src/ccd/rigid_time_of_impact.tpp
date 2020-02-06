#pragma once
#include "rigid_time_of_impact.hpp"

#include <ccd/interval_root_finder.hpp>
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
    int dim = bodyA.dim();
    assert(bodyB.dim() == dim);

    typedef Eigen::Matrix<Interval, Eigen::Dynamic, 1> VectorXI;

    const auto distance = [&](Interval t) -> Interval {
        physics::Pose<Interval> poseA = bodyA.pose.cast<Interval>()
            + displacementA.cast<Interval>() * toi;
        physics::Pose<Interval> poseB = bodyB.pose.cast<Interval>()
            + displacementB.cast<Interval>() * toi;

        VectorXI v0 = bodyA.world_vertex<Interval>(poseA, vertex_id);
        VectorXI v1
            = bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 0));
        VectorXI v2
            = bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 1));

        // TODO: Change this to define for 3D
        if (dim == 2) {
            VectorXI v2_v1 = v2 - v1;
            VectorXI normal(2);
            normal << v2_v1(1), v2_v1(0);
            normal.normalize();
            return (v0 - v1).transpose() * normal;
        } else {
            Interval v2_v1_sqrnorm = (v2 - v1).squaredNorm();
            // Ignored the case where 0 ∈ v2_v1_sqrnorm
            // if(zero_in(v2_v1_sqrnorm)) { ... }
            return ((v2 - v1).head<3>().cross((v1 - v0).head<3>()))
                       .squaredNorm()
                / v2_v1_sqrnorm;
        }
    };

    const auto is_point_along_edge = [&](Interval toi) {
        physics::Pose<Interval> poseA = bodyA.pose.cast<Interval>()
            + displacementA.cast<Interval>() * toi;
        physics::Pose<Interval> poseB = bodyB.pose.cast<Interval>()
            + displacementB.cast<Interval>() * toi;

        VectorXI v0 = bodyA.world_vertex<Interval>(poseA, vertex_id);
        VectorXI v1
            = bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 0));
        VectorXI v2
            = bodyB.world_vertex<Interval>(poseB, bodyB.edges(edge_id, 1));

        Interval alpha = (v1 - v0).dot((v2 - v1).normalized());
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

#include "rigid_time_of_impact.tpp"
