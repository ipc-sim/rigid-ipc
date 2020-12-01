// A spatial hash grid for rigid bodies with angular trajectories.
#pragma once

#include <ipc/spatial_hash/hash_grid.hpp>

#include <physics/rigid_body_assembler.hpp>

namespace ccd {

/// A spatial hash grid for rigid bodies with angular trajectories.
class RigidBodyHashGrid : public ipc::HashGrid {
public:
    void resize(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const std::vector<int>& body_ids,
        const double inflation_radius = 0.0);

    void addBodies(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const std::vector<int>& body_ids,
        const double inflation_radius = 0.0);

protected:
    void compute_vertices_intervals(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<Interval>& poses_t0,
        const physics::Poses<Interval>& poses_t1,
        // const std::vector<int>& body_ids,
        Eigen::MatrixXI& vertices,
        double inflation_radius = 0.0) const;

    int compute_vertices_intervals(
        const physics::RigidBody& body,
        const physics::Pose<Interval>& pose_t0,
        const physics::Pose<Interval>& pose_t1,
        // const std::vector<int>& body_ids,
        Eigen::MatrixXI& vertices,
        double inflation_radius = 0.0,
        const Interval& t = Interval(0, 1)) const;
};

} // namespace ccd
