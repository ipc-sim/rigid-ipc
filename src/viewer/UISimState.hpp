#pragma once

#include <memory> // shared_ptr

#include <spdlog/spdlog.h>

#include <igl/Timer.h>

#include <viewer/viewer_data.hpp>

// WARNING: Use an anonymous namespace when including gif.h to avoid duplicate
//          symbols.
namespace {
#include <gif.h>
}

#include "SimState.hpp"
#include <profiler.hpp>

namespace ipc::rigid {

class UISimState {
public:
    UISimState() = default;

    enum PlayerState { Playing = 0, Paused, TotalPlayerStatus };

    void init();
    void draw_menu();

    void launch(const std::string& inital_scene);
    void load_scene();
    void redraw_scene();
    void pre_draw_loop();
    void post_draw_loop();

    bool custom_key_pressed(unsigned int unicode_key, int modifiers);

    bool load(std::string scene_filename);
    bool save(std::string scene_filename)
    {
        return m_state.save_simulation(scene_filename);
    }

    bool save_obj_sequence(const std::string& dir_name)
    {
        bool success = m_state.save_obj_sequence(dir_name);
        m_state.problem_ptr->state(
            m_state.state_sequence[m_state.m_num_simulation_steps]);
        return success;
    }

    bool save_gltf(const std::string& filename)
    {
        return m_state.save_gltf(filename);
    }

    void start_recording(const std::string& filename);
    void end_recording();

    void reload();

    void simulation_step()
    {
        igl::Timer timer;
        timer.start();
        m_state.simulation_step();
        timer.stop();
        m_simulation_time += timer.getElapsedTime();
        m_state.save_simulation_step();
        redraw_scene();
    }

    void log_simulation_time()
    {
        spdlog::info("total_simulation_time={:g}s", m_simulation_time);
    }

    SimState m_state;
    PlayerState m_player_state = PlayerState::Paused;

    bool m_has_scene = false;            ///< @brief true if a scene is loaded
    bool m_bkp_had_collision = false;    ///< @brief breakpoint on collisions
    bool m_bkp_has_intersections = true; ///< @brief breakpoint on intersections
    bool m_bkp_optimization_failed = true; ///< @brief breakpoint on opt failure
    double m_interval_time = 0.0;          ///< @brief time within the interval
    bool m_show_vertex_data = false; ///< @brief show vertex data in the viewer

protected:
    void draw_io();
    void draw_simulation_player();
    void draw_settings();
    void draw_legends();

private:
    MeshData m_mesh_data;
    CoMData m_com_data;

    bool m_reloading_scene = false; ///< @brief is the scene is being reloaded?

    GifWriter m_gif_writer;
    uint32_t m_gif_delay = 1; //*10ms
    int m_gif_downscale = 1;
    bool m_is_gif_recording = false;
    bool m_scene_changed = false; ///< @brief scene changed since last draw

    double m_simulation_time = 0.0; ///< @brief total simulation time

    std::string inital_scene;
};

} // namespace ipc::rigid
