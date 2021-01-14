#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <physics/rigid_body.hpp>

namespace ccd {
namespace io {

    bool read_rb_scene_from_str(
        const std::string str, std::vector<physics::RigidBody>& rbs);

    bool read_rb_scene(
        const nlohmann::json& scene, std::vector<physics::RigidBody>& rbs);

} // namespace io
} // namespace ccd
