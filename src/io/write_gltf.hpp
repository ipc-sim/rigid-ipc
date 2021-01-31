#pragma once

#include <string>

#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

bool write_gltf(
    const std::string& filename,
    const RigidBodyAssembler& bodies,
    const std::vector<PosesD>& poses,
    double timestep,
    bool embed_buffers = true,
    bool write_binary = true,
    bool prettyPrint = true);

} // namespace ipc::rigid
