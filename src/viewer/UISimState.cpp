#include "UISimState.hpp"
#include <logger.hpp>

namespace ccd {

UISimState::UISimState()
    : m_player_state(PlayerState::Paused)
    , m_has_scene(false)
    , m_bkp_had_collision(false)
    , m_bkp_has_collision(true)
    , m_log_level(spdlog::level::info)
    , m_interval_time(0.0)
    , m_show_as_delta(true)
    , m_show_next_step(true)
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
    collision_force_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer,
        Eigen::RowVector3d(236, 240, 241) / 255.0); // #ecf0f1 - CLOUDS WHITE
    collision_force_data->data().show_overlay = false;

    viewer->append_mesh();
    gradient_data = std::make_unique<igl::opengl::VectorFieldData>(_viewer,
        Eigen::RowVector3d(243, 156, 18) / 255.0); // #f39c12 - ORANGE

    viewer->append_mesh();
    edges_data = std::make_unique<igl::opengl::GraphData>(_viewer,
        Eigen::RowVector3d(231.0, 76, 60) / 255.0); // #e74c3c - ALIZARIN RED

    viewer->append_mesh();
    next_displacement_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer,
        Eigen::RowVector3d(46, 204, 113) / 255.0); // #2ecc71 - EMERALD GREEN

    viewer->append_mesh();
    initial_displacement_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer, Eigen::RowVector3d(46, 204, 113) / 255.0); // #2ecc71

    viewer->append_mesh();
    final_displacement_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer, Eigen::RowVector3d(155, 89, 182) / 255.0); // #9b59b6 AMETHIST

    viewer->append_mesh();
    velocity_data = std::make_unique<igl::opengl::VectorFieldData>(_viewer,
        Eigen::RowVector3d(241, 196, 15) / 255.0); // #f1c40f SUN FLOWER
    velocity_data->data().show_overlay = false;

    // keep this last!
    viewer->append_mesh();
    grid_data = std::make_unique<igl::opengl::ScalarFieldData>(_viewer,
        /*color_inf=*/Eigen::RowVector3d::Zero(),
        Eigen::RowVector3d(0.3, 0.3, 0.5)); // igl default
    grid_data->data().show_faces = false;

    datas_.emplace("edges", edges_data);
    datas_.emplace("gradient", gradient_data);
    datas_.emplace("next displ.", next_displacement_data);
    datas_.emplace("initial displ.", initial_displacement_data);
    datas_.emplace("final displ.", final_displacement_data);
    datas_.emplace("velocity", velocity_data);
    datas_.emplace("F_c", collision_force_data);
    datas_.emplace("grid", grid_data);

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
    auto q0 = m_state.problem_ptr->vertices_prev();
    auto q1 = m_state.problem_ptr->vertices();
    auto q1_coll = m_state.problem_ptr->vertices_collision();
    auto q2 = m_state.problem_ptr->vertices_next(m_state.m_timestep_size);
    auto v1 = m_state.problem_ptr->velocities(
        m_show_as_delta, m_state.m_timestep_size);
    auto fc = m_state.problem_ptr->collision_force(
        m_show_as_delta, m_state.m_timestep_size);

    edges_data->data().show_vertid = false;
    edges_data->set_graph(q1, m_state.problem_ptr->edges());
    edges_data->set_vertex_data(m_state.problem_ptr->particle_dof_fixed());
    edges_data->data().point_size = 10 * pixel_ratio();

    next_displacement_data->set_vector_field(q1, q2 - q1);
    initial_displacement_data->set_vector_field(q0, q1_coll - q0);
    final_displacement_data->set_vector_field(q0, q1 - q0);

    velocity_data->set_vector_field(q1, v1);
    collision_force_data->set_vector_field(q1, -fc);

    Eigen::MatrixXd grad_f;
    auto x1 = m_state.problem_ptr->dof_positions();
    m_state.get_collision_gradient(grad_f);
    gradient_data->set_vector_field(x1, -grad_f);

    grid_data->set_mesh(m_state.grid_V, m_state.grid_F);

    viewer->core().align_camera_center(edges_data->mV, edges_data->mE);
    m_has_scene = true;
    m_player_state = PlayerState::Paused;
    m_interval_time = 0.0;
}

void UISimState::redraw_scene()
{
    auto q0 = m_state.problem_ptr->vertices_prev();
    auto q1 = m_state.problem_ptr->vertices();
    auto q1_coll = m_state.problem_ptr->vertices_collision();
    auto q2 = m_state.problem_ptr->vertices_next(m_state.m_timestep_size);
    auto v1 = m_state.problem_ptr->velocities(
        m_show_as_delta, m_state.m_timestep_size);
    auto fc = m_state.problem_ptr->collision_force(
        m_show_as_delta, m_state.m_timestep_size);

    if (m_show_next_step) {
        edges_data->update_graph(q1 + m_interval_time * (q2 - q1));
    } else {
        edges_data->update_graph(q0 + m_interval_time * (q1 - q0));
    }

    next_displacement_data->update_vector_field(q1, q2 - q1);
    initial_displacement_data->update_vector_field(q0, q1_coll - q0);
    final_displacement_data->update_vector_field(q0, q1 - q0);

    Eigen::MatrixXd grad_f;
    auto x1 = m_state.problem_ptr->dof_positions();
    m_state.get_collision_gradient(grad_f);
    gradient_data->update_vector_field(x1, -grad_f);

    velocity_data->update_vector_field(q1, v1);
    collision_force_data->update_vector_field(q1, -fc);
    redraw_scalar_fields();
}

void UISimState::redraw_scalar_fields()
{
    // Computing this is expensive so we do it on demand
    if (grid_data->visibility()) {
        Eigen::VectorXd fx;
        m_state.get_collision_functional_field(fx);
        grid_data->set_vertex_data(fx);
    }
}

bool UISimState::pre_draw_loop()
{
    if (m_player_state == PlayerState::Playing) {
        simulation_step();
        bool breakpoint = m_bkp_had_collision && m_state.m_step_had_collision;
        breakpoint = breakpoint
            || (m_bkp_has_collision && m_state.m_step_has_collision);

        if (breakpoint) {
            m_player_state = PlayerState::Paused;
        }
    }
    return false;
}
} // namespace ccd
