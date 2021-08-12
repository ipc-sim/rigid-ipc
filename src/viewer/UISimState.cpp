#include "UISimState.hpp"
#include <logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>

#include <physics/rigid_body_problem.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

UISimState::UISimState()
    : m_player_state(PlayerState::Paused)
    , m_has_scene(false)
    , m_bkp_had_collision(false)
    , m_bkp_has_intersections(true)
    , m_bkp_optimization_failed(true)
    , m_interval_time(0.0)
    , m_show_vertex_data(false)
    , m_reloading_scene(false)
    , m_scene_changed(false)
    , m_simulation_time(0)
{
}

void UISimState::launch(const std::string& inital_scene)
{
    m_viewer.plugins.push_back(this);
    m_viewer.core().set_rotation_type(
        igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    m_viewer.core().orthographic = true;
    m_viewer.core().is_animating = true;
    m_viewer.core().lighting_factor = 0.0;
    m_viewer.core().animation_max_fps = 120.0;
    this->inital_scene = inital_scene;
    m_viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) {
        return pre_draw_loop();
    };
    m_viewer.callback_post_draw = [&](igl::opengl::glfw::Viewer&) {
        return post_draw_loop();
    };
    m_viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&,
                                        unsigned int unicode_key,
                                        int modifiers) {
        return custom_key_pressed(unicode_key, modifiers);
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
    mesh_data->visibility(true);

    viewer->append_mesh();
    velocity_data = std::make_unique<igl::opengl::VectorFieldData>(
        _viewer,
        Eigen::RowVector3d(0xF1, 0xC4, 0x0F) / 0xFF); // #F1C40F - SUN FLOWER
    velocity_data->visibility(false);

    viewer->append_mesh();
    com_data = std::make_unique<igl::opengl::CoMData>(_viewer);
    com_data->visibility(true);

    datas_.emplace_back("edges", mesh_data);
    datas_.emplace_back("velocity", velocity_data);
    datas_.emplace_back("body-frame", com_data);

    load(this->inital_scene);
}

void UISimState::load_scene()
{
    Eigen::MatrixXd q = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v =
        m_state.problem_ptr->velocities() * m_state.problem_ptr->timestep();

    // mesh_data->data().show_vertid = false;
    mesh_data->set_mesh(
        q, m_state.problem_ptr->edges(), m_state.problem_ptr->faces());
    Eigen::VectorXi vertex_type;
    if (m_state.problem_ptr->is_rb_problem()) {
        const auto& bodies =
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler;
        vertex_type.resize(bodies.num_vertices());
        int start_i = 0;
        for (const auto& body : bodies.m_rbs) {
            vertex_type.segment(start_i, body.vertices.rows())
                .setConstant(int(body.type));
            start_i += body.vertices.rows();
        }
    } else {
        vertex_type =
            m_state.problem_ptr->vertex_dof_fixed().rowwise().all().cast<int>();
        vertex_type *= 2; // 0 is static, 1 is kinematic, 2 is dynamic
    }
    mesh_data->set_vertex_data(
        m_state.problem_ptr->vertex_dof_fixed(), vertex_type);
    mesh_data->data().point_size = 0 * pixel_ratio();

    velocity_data->set_vector_field(q, v);

    if (m_state.problem_ptr->is_rb_problem()) {
        com_data->set_coms(
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler.rb_poses());
    }

    m_has_scene = true;
    m_player_state = PlayerState::Paused;
    m_interval_time = 0.0;

    m_simulation_time = 0;

    // Do not change the view setting upon reload
    if (m_reloading_scene) {
        return;
    }

    int dim = q.cols();
    m_viewer.core().trackball_angle.setIdentity();
    m_viewer.core().set_rotation_type(
        dim == 2 ? igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION
                 : igl::opengl::ViewerCore::
                       ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);
    m_viewer.core().orthographic = dim == 2;
    m_viewer.core().lighting_factor = 0.0; // dim == 2 ? 0.0 : 1.0;
    // mesh_data->data().set_face_based(true);
    m_viewer.core().align_camera_center(
        q, dim == 2.0 ? mesh_data->mE : mesh_data->mF);

    // Default colors
    // background_color << 0.3f, 0.3f, 0.5f, 1.0f;

    // Camera parameters
    m_viewer.core().camera_zoom = 1.0f;
    m_viewer.core().camera_translation << 0, 0, 0;
}

void UISimState::redraw_scene()
{
    Eigen::MatrixXd q1 = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v1 =
        m_state.problem_ptr->velocities() * m_state.problem_ptr->timestep();

    mesh_data->update_vertices(q1);
    velocity_data->update_vector_field(q1, v1);
    if (m_state.problem_ptr->is_rb_problem()) {
        const auto& bodies =
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler;
        Eigen::VectorXi vertex_type(bodies.num_vertices());
        int start_i = 0;
        for (const auto& body : bodies.m_rbs) {
            vertex_type.segment(start_i, body.vertices.rows())
                .setConstant(int(body.type));
            start_i += body.vertices.rows();
        }
        mesh_data->set_vertex_data(
            m_state.problem_ptr->vertex_dof_fixed(), vertex_type);
    }

    if (m_state.problem_ptr->is_rb_problem()) {
        com_data->set_coms(
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler.rb_poses());
    }
}

bool UISimState::pre_draw_loop()
{
    size_t last_save_state = m_state.state_sequence.size() - 1;
    if (m_state.m_num_simulation_steps > last_save_state) {
        m_state.m_num_simulation_steps = last_save_state;
        m_player_state = PlayerState::Paused;
        replaying = false;
    }
    if (replaying) {
        m_state.problem_ptr->state(
            m_state.state_sequence[m_state.m_num_simulation_steps]);
        redraw_scene();
        m_scene_changed = true;
        if (m_player_state == PlayerState::Playing) {
            m_state.m_num_simulation_steps++;
        } else if (m_state.m_num_simulation_steps >= last_save_state) {
            m_state.m_num_simulation_steps = last_save_state;
            m_player_state = PlayerState::Paused;
            replaying = false;
        }
    } else if (m_player_state == PlayerState::Playing) {
        simulation_step();

        bool breakpoint = (m_bkp_had_collision && m_state.m_step_had_collision)
            || (m_bkp_has_intersections && m_state.m_step_has_intersections)
            || (m_bkp_optimization_failed
                && !m_state.problem_ptr->opt_result.success)
            || (m_state.m_max_simulation_steps >= 0
                && m_state.m_num_simulation_steps
                    >= m_state.m_max_simulation_steps);
        if (breakpoint) {
            m_player_state = PlayerState::Paused;
            log_simulation_time();
        }
        m_scene_changed = true;
    }
    return false;
}

void UISimState::save_screenshot(const std::string& filename)
{
    if (filename == "") {
        return;
    }

    int width, height;
    get_window_dimensions(width, height);
    using MatrixXuc = MatrixX<unsigned char>;
    // Allocate temporary buffers for image
    MatrixXuc R(width, height), G(width, height), B(width, height),
        A(width, height);

    // Draw the scene in the buffers
    // m_viewer.core().draw_buffer(mesh_data->data(), false, R, G, B, A);
    // m_viewer.core().draw_buffer(velocity_data->data(), true, R, G, B, A);
    // m_viewer.core().draw_buffer(com_data->data(), true, R, G, B, A);

    std::vector<unsigned char> data(4 * width * height);
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data.data());
    // img->flip();
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            R(i, j) = data[4 * (i + j * width) + 0];
            G(i, j) = data[4 * (i + j * width) + 1];
            B(i, j) = data[4 * (i + j * width) + 2];
            A(i, j) = data[4 * (i + j * width) + 3];
        }
    }
    bool success = igl::png::writePNG(R, G, B, A, filename);

    if (!success) {
        spdlog::error("Unable to save screenshot to {}", filename);
    }
}

void UISimState::start_recording(const std::string& filename)
{
    if (filename == "") {
        return;
    }

    int width, height;
    get_window_dimensions(width, height);
    width = static_cast<int>(m_gif_scale * width);
    height = static_cast<int>(m_gif_scale * height);
    GifBegin(&m_gif_writer, filename.c_str(), width, height, m_gif_delay);
    m_scene_changed = true;
    // post_draw_loop();
    m_is_gif_recording = true;
}

bool UISimState::post_draw_loop()
{
    if (m_scene_changed && m_is_gif_recording) {
        int width, height;
        get_window_dimensions(width, height);
        width = static_cast<int>(m_gif_scale * width);
        height = static_cast<int>(m_gif_scale * height);

        using MatrixXuc = MatrixX<unsigned char>;
        // Allocate temporary buffers for image
        MatrixXuc R(width, height), G(width, height), B(width, height),
            A(width, height);

        // Draw the scene in the buffers
        m_viewer.core().draw_buffer(mesh_data->data(), false, R, G, B, A);

        std::vector<uint8_t> img(width * height * 4);
        for (int rowI = 0; rowI < width; rowI++) {
            for (int colI = 0; colI < height; colI++) {
                int indStart = (rowI + (height - 1 - colI) * width) * 4;
                img[indStart] = R(rowI, colI);
                img[indStart + 1] = G(rowI, colI);
                img[indStart + 2] = B(rowI, colI);
                img[indStart + 3] = A(rowI, colI);
            }
        }
        GifWriteFrame(&m_gif_writer, img.data(), width, height, m_gif_delay);
    }
    m_scene_changed = false;
    return false;
}

void UISimState::end_recording()
{
    GifEnd(&m_gif_writer);
    m_is_gif_recording = false;
}

bool UISimState::custom_key_pressed(unsigned int unicode_key, int modifiers)
{
    if (mesh_data == nullptr) {
        return false;
    }

    switch (unicode_key) {
    case 'F':
    case 'f': {
        // mesh_data->data().set_face_based(!mesh_data->data().face_based);
        return true;
    }
    case 'L':
    case 'l': {
        bool toggle = mesh_data->visibility();
        mesh_data->visibility(!toggle);
        return true;
    }
    case 'T':
    case 't': {
        m_viewer.core().toggle(mesh_data->data().show_faces);
        return true;
    }
    case ' ': {
        m_player_state = m_player_state == PlayerState::Playing
            ? PlayerState::Paused
            : PlayerState::Playing;
        log_simulation_time();
        return true;
    }
    default:
        break; // do nothing
    }
    return false;
}

} // namespace ipc::rigid
