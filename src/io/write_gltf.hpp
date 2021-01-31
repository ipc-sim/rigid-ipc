#pragma once

#include <string>

#include <physics/rigid_body_assembler.hpp>

namespace ccd {
namespace io {

    bool write_gltf(
        const std::string& filename,
        const physics::RigidBodyAssembler& bodies,
        const std::vector<physics::Poses<double>>& poses,
        double timestep,
        bool embed_buffers = true,
        bool write_binary = true,
        bool prettyPrint = true);

}
} // namespace ccd
