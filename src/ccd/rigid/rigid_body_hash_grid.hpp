// A spatial hash grid for rigid bodies with angular trajectories.
#pragma once

#include <ipc/broad_phase/hash_grid.hpp>

#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

/// A spatial hash grid for rigid bodies with angular trajectories.
class RigidBodyHashGrid : public HashGrid {
public:
    /// Resize to fit a static scene
    void resize(
        const RigidBodyAssembler& bodies,
        const PosesD& poses,
        const std::vector<std::pair<int, int>>& body_pairs,
        const double inflation_radius = 0.0);

    /// Resize to fit a dynamic scene
    void resize(
        const RigidBodyAssembler& bodies,
        const PosesD& poses_t0,
        const PosesD& poses_t1,
        const std::vector<std::pair<int, int>>& body_pairs,
        const double inflation_radius = 0.0);

    /// Add static bodies
    void addBodies(
        const RigidBodyAssembler& bodies,
        const PosesD& poses,
        const std::vector<std::pair<int, int>>& body_pairs,
        const double inflation_radius = 0.0);

    /// Add dynamic bodies
    void addBodies(
        const RigidBodyAssembler& bodies,
        const PosesD& poses_t0,
        const PosesD& poses_t1,
        const std::vector<std::pair<int, int>>& body_pairs,
        const double inflation_radius = 0.0);

protected:
    void compute_vertices_intervals(
        const RigidBodyAssembler& bodies,
        const Poses<Interval>& poses_t0,
        const Poses<Interval>& poses_t1,
        const std::vector<int>& body_ids,
        MatrixXI& vertices,
        double inflation_radius = 0.0) const;

    int compute_vertices_intervals(
        const RigidBody& body,
        const Pose<Interval>& pose_t0,
        const Pose<Interval>& pose_t1,
        MatrixXI& vertices,
        double inflation_radius = 0.0,
        const Interval& t = Interval(0, 1),
        int force_subdivision = 0) const;
};

} // namespace ipc::rigid
