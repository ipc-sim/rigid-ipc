#include "UISimState.hpp"
#include <logger.hpp>

#include <spdlog/sinks/stdout_color_sinks.h>

#include <physics/rigid_body_problem.hpp>
#include <utils/eigen_ext.hpp>

#include <polyscope/polyscope.h>

namespace ps = polyscope;

namespace ipc::rigid {

void UISimState::launch(const std::string& inital_scene)
{
    this->inital_scene = inital_scene;

    ps::options::programName = "Rigid IPC";
    ps::options::giveFocusOnShow = true;

    ps::init();
    ps::options::groundPlaneMode = ps::GroundPlaneMode::None;

    this->init();

    ps::state::userCallback = [this]() {
        draw_menu();
        pre_draw_loop();
        post_draw_loop();
    };

    ps::show();
}

void UISimState::init() { load(this->inital_scene); }

bool UISimState::load(std::string scene_filename)
{
    if (scene_filename != "" && m_state.load_scene(scene_filename)) {
        load_scene();
        return true;
    }
    return false;
}

void UISimState::reload()
{
    m_reloading_scene = true;
    m_state.reload_scene();
    load_scene();
    m_reloading_scene = false;
}

void UISimState::load_scene()
{
    Eigen::MatrixXd q = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v = m_state.problem_ptr->velocities();

    m_mesh_data.set_mesh(
        q, m_state.problem_ptr->edges(), m_state.problem_ptr->faces(),
        m_state.problem_ptr->codim_edges_to_edges(),
        m_state.problem_ptr->codim_vertices_to_vertices());
    m_mesh_data.update_velocities(v);

    Eigen::VectorXi vertex_types, vertex_ids;
    if (m_state.problem_ptr->is_rb_problem()) {
        const auto& bodies =
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler;
        vertex_types.resize(bodies.num_vertices());
        vertex_ids.resize(bodies.num_vertices());
        int start_i = 0;
        int i = 0;
        for (const auto& body : bodies.m_rbs) {
            vertex_types.segment(start_i, body.vertices.rows())
                .setConstant(int(body.type));
            vertex_ids.segment(start_i, body.vertices.rows()).setConstant(i++);
            start_i += body.vertices.rows();
        }

    } else {
        vertex_types =
            m_state.problem_ptr->vertex_dof_fixed().rowwise().all().cast<int>();
        vertex_types *= 2; // 0 is static, 1 is kinematic, 2 is dynamic
        vertex_ids = m_state.problem_ptr->group_ids();
    }
    m_mesh_data.set_vertex_types(vertex_types);
    m_mesh_data.set_vertex_ids(vertex_ids);

    if (m_state.problem_ptr->is_rb_problem()) {
        m_com_data.set_coms(
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
    // m_viewer.core().trackball_angle.setIdentity();
    if (dim == 2) {
        ps::view::setNavigateStyle(ps::NavigateStyle::Planar);
        ps::view::projectionMode = ps::ProjectionMode::Orthographic;
    } else {
        ps::view::setNavigateStyle(ps::NavigateStyle::Turntable);
        ps::view::projectionMode = ps::ProjectionMode::Perspective;
    }
    // m_viewer.core().align_camera_center(
    //     q, dim == 2.0 ? mesh_data->mE : mesh_data->mF);

    // Default colors
    // if (q.cols() == 2) {
    ps::view::bgColor = { { 0.0f, 0.0f, 0.0f, 0.0f } }; // Black for 2D
    // } else {
    // ps::view::bgColor = { { 1.0f, 1.0f, 1.0f, 0.0f } }; // White for 3D
    // }
    // ps::view::bgColor = { { 0.3f, 0.3f, 0.5f, 0.0f } };

    // Camera parameters
    // m_viewer.core().camera_zoom = 1.0f;
    // m_viewer.core().camera_translation << 0, 0, 0;
}

void UISimState::pre_draw_loop()
{
    size_t last_save_state = m_state.state_sequence.size() - 1;
    if (m_state.m_num_simulation_steps > last_save_state) {
        m_state.m_num_simulation_steps = last_save_state;
        m_player_state = PlayerState::Paused;
    }

    if (m_player_state == PlayerState::Playing) {
        if (m_state.m_num_simulation_steps < last_save_state) {
            // Replaying the simulation
            m_state.problem_ptr->state(
                m_state.state_sequence[++m_state.m_num_simulation_steps]);
            redraw_scene();
            if (m_state.m_num_simulation_steps >= last_save_state) {
                m_state.m_num_simulation_steps = last_save_state;
                m_player_state = PlayerState::Paused;
            }
        } else {
            simulation_step(); // calls redraw_scene()

            bool breakpoint =
                (m_bkp_had_collision && m_state.m_step_had_collision)
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
        }
        m_scene_changed = true;
    }
}

void UISimState::redraw_scene()
{
    Eigen::MatrixXd q1 = m_state.problem_ptr->vertices();
    Eigen::MatrixXd v1 = m_state.problem_ptr->velocities();

    m_mesh_data.update_vertices(q1);
    m_mesh_data.update_velocities(v1);

    if (m_state.problem_ptr->is_rb_problem()) {
        m_com_data.update_coms(
            std::dynamic_pointer_cast<RigidBodyProblem>(m_state.problem_ptr)
                ->m_assembler.rb_poses());
    }
}

void UISimState::post_draw_loop()
{
    if (m_scene_changed && m_is_gif_recording) {
        int width = ps::view::bufferWidth, height = ps::view::bufferHeight;
        width = width >> m_gif_downscale;
        height = height >> m_gif_downscale;

        const std::vector<unsigned char> buffer = ps::screenshotToBuffer(
            /*transparentBG=*/true);

        // Flip the buffer vertically and downscale it
        std::vector<uint8_t> img(width * height * 4);
        const int delta = 1 << m_gif_downscale;
        for (int x = 0; x < ps::view::bufferWidth; x += delta) {
            for (int y = 0; y < ps::view::bufferHeight; y += delta) {
                // int buffer_i = (y * ps::view::bufferWidth + x) * 4;
                int img_i = ((height - 1 - (y >> m_gif_downscale)) * width
                             + (x >> m_gif_downscale))
                    * 4;
                for (int d = 0; d < 4; ++d) {
                    img[img_i + d] = 0; // Initialize to zero
                    // Average the pixels in the downscaled region
                    for (int i = 0; i < delta; ++i) {
                        for (int j = 0; j < delta; ++j) {
                            int buffer_i =
                                ((y + i) * ps::view::bufferWidth + (x + j)) * 4;
                            img[img_i + d] +=
                                buffer[buffer_i + d] / (delta * delta);
                        }
                    }
                }
            }
        }

        GifWriteFrame(&m_gif_writer, img.data(), width, height, m_gif_delay);
    }
    m_scene_changed = false;
}

void UISimState::start_recording(const std::string& filename)
{
    if (filename == "") {
        return;
    }

    int width = ps::view::bufferWidth, height = ps::view::bufferHeight;
    width = width >> m_gif_downscale;
    height = height >> m_gif_downscale;
    GifBegin(&m_gif_writer, filename.c_str(), width, height, m_gif_delay);
    m_is_gif_recording = true;
}

void UISimState::end_recording()
{
    GifEnd(&m_gif_writer);
    m_is_gif_recording = false;
}

} // namespace ipc::rigid
