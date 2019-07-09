#include "rigid_body_view.hpp"

#include <imgui/imgui.h>
#include <viewer/imgui_ext.hpp>

namespace ccd {

bool rigid_body_menu(physics::RigidBody& rb)
{
    bool update = false;
    update = ImGui::DragDouble("x pos.##rigid-bodies", &rb.position(0), 0.1)
        || update;
    update = ImGui::DragDouble("y pos.##rigid-bodies", &rb.position(1), 0.1)
        || update;
    double theta_deg = rb.position(2) * 180 / M_PI;
    if (ImGui::DragDouble("theta(deg)##rigid-bodies", &theta_deg, 0.1)) {
        update = true;
        rb.position(2) = theta_deg * M_PI / 180;
    }

    update
        = ImGui::InputDouble("x vel.##rigid-bodies", &rb.velocity(0)) || update;
    update
        = ImGui::InputDouble("y vel.##rigid-bodies", &rb.velocity(1)) || update;
    double omega_deg = rb.velocity(2) * 180 / M_PI;
    if (ImGui::DragDouble("omega(deg)##rigid-bodies", &omega_deg, 0.1)) {
        update = true;
        rb.velocity(2) = omega_deg * M_PI / 180;
    }

    return update;
}

bool rigid_body_system_menu(
    physics::RigidBodySystem& rbs, const std::vector<int>& selected_points)
{
    if (selected_points.size() > 0) {
        size_t selected = rbs.vertex_to_body_map(selected_points[0]);
        return rigid_body_menu(rbs.rigid_bodies[selected]);
    }
    return false;
}

} // namespace ccd
