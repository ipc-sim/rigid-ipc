#pragma once

#include <ipc/spatial_hash/collision_candidate.hpp>
#include <physics/pose.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {

void save_ccd_candidate(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeVertexCandidate& ev_candidate);

void save_ccd_candidate(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::FaceVertexCandidate& fv_candidate);

void save_ccd_candidate(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const ipc::EdgeEdgeCandidate& ee_candidate);

} // namespace ccd
