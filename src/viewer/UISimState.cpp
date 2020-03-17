#include "UISimState.hpp"
#include <logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>

namespace ccd {

UISimState::UISimState()
    : m_player_state(PlayerState::Paused)
    , m_has_scene(false)
    , m_bkp_had_collision(false)
    , m_bkp_has_intersections(true)
    , m_log_level(spdlog::level::info)
    , m_interval_time(0.0)
    , m_show_vertex_data(false)
    , m_reloading_scene(false)
{
}

void UISimState::launch()
{
    m_viewer.plugins.push_back(this);
    m_viewer.core().set_rotation_type(
        igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    m_viewer.core().orthographic = true;
    m_viewer.core().is_animating = true;
    m_viewer.core().lighting_factor = 0.0;
    m_viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) {
        return pre_draw_loop();
    };
    m_viewer.launch(
        /*resizable=*/true, /*fullscreen=*/false,
        /*name=*/"Simulation");
}

void UISimState::init(igl::opengl::glfw::Viewer* _viewer)
{
    Super::init(_viewer);
    viewer->data().clear();

    viewer->append_mesh();
    mesh_data = std::make_unique<igl::opengl::MeshData>(
        _viewer,
        Eigen::RowVector3d(0xE7, 0x4C, 0x3C) / 0xFF); // #E74C3C - ALIZARIN RED
    viewer->append_mesh();
    velocity_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer,
        Eigen::RowVector3d(0xF1, 0xC4, 0x0F) / 0xFF); // #F1C40F SUN FLOWER
    velocity_data->data().show_overlay = false;

    datas_.emplace("edges", mesh_data);
    datas_.emplace("velocity", velocity_data);

    for (auto it = datas_.begin(); it != datas_.end(); ++it) {
        data_names_.push_back(it->first);
    }
}

std::shared_ptr<igl::opengl::ViewerDataExt>
UISimState::get_data(const std::string& dataname) const
{
    auto it = datas_.find(dataname);
    assert(it != datas_.end());
    return it->second;
}

void UISimState::load_scene()
{
    Eigen::MatrixXd q = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v =
        m_state.problem_ptr->velocities() * m_state.m_timestep_size;

    mesh_data->data().show_vertid = false;
    mesh_data->set_mesh(
        q, m_state.problem_ptr->edges(), m_state.problem_ptr->faces());
    mesh_data->set_vertex_data(m_state.problem_ptr->particle_dof_fixed());
    mesh_data->data().point_size = 0 * pixel_ratio();

    velocity_data->set_vector_field(q, v);

    m_has_scene = true;
    m_player_state = PlayerState::Paused;
    m_interval_time = 0.0;

    // Do not change the view setting upon reload
    if (m_reloading_scene) {
        return;
    }

    int dim = q.cols();
    m_viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
    m_viewer.core().set_rotation_type(
        dim == 2 ? igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION
                 : igl::opengl::ViewerCore::
                       ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);
    m_viewer.core().orthographic = dim == 2;
    m_viewer.core().lighting_factor = float(dim == 2);
    viewer->core().align_camera_center(mesh_data->mV, mesh_data->mE);

    // Default colors
    // background_color << 0.3f, 0.3f, 0.5f, 1.0f;

    // Camera parameters
    if (dim == 3) {
        m_viewer.core().camera_base_zoom = 1.0f;
        m_viewer.core().camera_zoom = 1.0f;
        m_viewer.core().camera_view_angle = 45.0;
        m_viewer.core().camera_base_translation << 0, 0, 0;
        m_viewer.core().camera_translation << 0, 0, 0;
        m_viewer.core().camera_eye << 0, 0, 5;
        m_viewer.core().camera_center << 0, 0, 0;
        m_viewer.core().camera_up << 0, 1, 0;
    }
}

void UISimState::redraw_scene()
{
    Eigen::MatrixXd q1 = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v1 =
        m_state.problem_ptr->velocities() * m_state.m_timestep_size;

    mesh_data->update_vertices(q1);
    velocity_data->update_vector_field(q1, v1);
}

bool UISimState::pre_draw_loop()
{
    if (m_player_state == PlayerState::Playing) {
        simulation_step();
        bool breakpoint = m_bkp_had_collision && m_state.m_step_had_collision;
        breakpoint = breakpoint
            || (m_bkp_has_intersections && m_state.m_step_has_collision);
        if (m_state.m_max_simulation_steps > -1) {
            breakpoint = breakpoint
                || (m_state.m_num_simulation_steps
                    >= m_state.m_max_simulation_steps);
        }
        if (breakpoint) {
            m_player_state = PlayerState::Paused;
        }
    }
    return false;
}
} // namespace ccd
