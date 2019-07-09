#include "UISimState.hpp"
#include <logger.hpp>

namespace ccd {

UISimState::UISimState()
    : m_player_state(PlayerState::Paused)
    , m_has_scene(false)
    , m_bkp_had_collision(false)
    , m_bkp_has_collision(true)
    , m_log_level(spdlog::level::info)
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
    edges_data = std::make_unique<igl::opengl::ViewerDataExt>(viewer->data());
}

void UISimState::load_scene()
{
    edges_data->m_data.show_vertid = true;
    edges_data->set_graph(m_state.m_problem.m_assembler.world_vertices(),
        m_state.m_problem.m_assembler.m_edges,
        Eigen::RowVector3d(1.0, 0.0, 0.0)); // #ff0000
    edges_data->m_data.point_size = 10 * pixel_ratio();

    viewer->core().align_camera_center(edges_data->mV, edges_data->mE);
    m_has_scene = true;
}

void UISimState::redraw_scene()
{
    edges_data->update_graph(m_state.m_problem.m_assembler.world_vertices());
}

bool UISimState::pre_draw_loop()
{
    if (m_player_state == PlayerState::Playing) {
        simulation_step();
        bool breakpoint = m_bkp_had_collision && m_state.m_step_had_collision;
        breakpoint = breakpoint || (m_bkp_has_collision && m_state.m_step_has_collision);

        if (breakpoint){
            m_player_state = PlayerState::Paused;
        }
    }
    return false;
}
} // namespace ccd
