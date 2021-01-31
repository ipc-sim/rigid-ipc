#pragma once

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include <physics/rigid_body.hpp>

namespace ipc::rigid {

bool read_rb_scene_from_str(
    const std::string str, std::vector<RigidBody>& rbs);

bool read_rb_scene(
    const nlohmann::json& scene, std::vector<RigidBody>& rbs);

} // namespace ipc::rigid
