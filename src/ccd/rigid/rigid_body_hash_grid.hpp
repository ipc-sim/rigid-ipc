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
};

} // namespace ccd
