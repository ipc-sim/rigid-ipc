#include "UISimState.hpp"
#include <logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>

namespace ccd {

UISimState::UISimState()
    : m_player_state(PlayerState::Paused)
    , m_has_scene(false)
    , m_bkp_had_collision(false)
    , m_bkp_has_collision(true)
    , m_log_level(spdlog::level::info)
    , m_interval_time(0.0)
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
    m_viewer.callback_pre_draw
        = [&](igl::opengl::glfw::Viewer&) { return pre_draw_loop(); };
    m_viewer.launch();
}

void UISimState::init(igl::opengl::glfw::Viewer* _viewer)
{

    Super::init(_viewer);
    viewer->data().clear();

    viewer->append_mesh();
    edges_data = std::make_unique<igl::opengl::GraphData>(_viewer,
        Eigen::RowVector3d(231.0, 76, 60) / 255.0); // #e74c3c - ALIZARIN RED    
    datas_.emplace("edges", edges_data);

    for (auto it = datas_.begin(); it != datas_.end(); ++it) {
        data_names_.push_back(it->first);
    }
}

std::shared_ptr<igl::opengl::ViewerDataExt> UISimState::get_data(
    const std::string& dataname) const
{
    auto it = datas_.find(dataname);
    assert(it != datas_.end());
    return it->second;
}

void UISimState::load_scene()
{

    auto q = m_state.problem_ptr->vertices();

    edges_data->data().show_vertid = false;
    edges_data->set_graph(q, m_state.problem_ptr->edges());
    edges_data->set_vertex_data(m_state.problem_ptr->particle_dof_fixed());
    edges_data->data().point_size = 10 * pixel_ratio();

    viewer->core().align_camera_center(edges_data->mV, edges_data->mE);
    m_has_scene = true;
    m_player_state = PlayerState::Paused;
    m_interval_time = 0.0;
}

void UISimState::redraw_scene()
{
    auto q1 = m_state.problem_ptr->vertices();
    edges_data->update_graph(q1);
}

bool UISimState::pre_draw_loop()
{
    if (m_player_state == PlayerState::Playing) {
        simulation_step();
        bool breakpoint = m_bkp_had_collision && m_state.m_step_had_collision;
        breakpoint = breakpoint
            || (m_bkp_has_collision && m_state.m_step_has_collision);
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
