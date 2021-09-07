#pragma once

#include <ipc/broad_phase/collision_candidate.hpp>
#include <physics/pose.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& ev_candidate);

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& fv_candidate);

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& ee_candidate);

} // namespace ipc::rigid
