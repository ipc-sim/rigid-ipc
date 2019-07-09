#pragma once

#include <physics/rigid_body.hpp>
#include <physics/rigid_body_system.hpp>
namespace ccd {

    bool rigid_body_menu(physics::RigidBody& rb);
    bool rigid_body_system_menu(physics::RigidBodySystem& rbs, const std::vector<int>& selected_points);
}
