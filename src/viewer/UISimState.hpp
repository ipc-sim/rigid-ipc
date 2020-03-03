#pragma once

#include <memory> // shared_ptr

#include <spdlog/spdlog.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <viewer/igl_viewer_ext.hpp>

#include "SimState.hpp"
#include <profiler.hpp>

namespace ccd {

class UISimState : public igl::opengl::glfw::imgui::ImGuiMenu {
    typedef igl::opengl::glfw::imgui::ImGuiMenu Super;

public:
    UISimState();
    ~UISimState() override {}

    enum PlayerState { Playing = 0, Paused, TotalPlayerStatus };

    virtual void init(igl::opengl::glfw::Viewer* _viewer) override;
    virtual void draw_menu() override;
    // virtual bool mouse_down(int button, int modifier) override;
    // virtual bool key_pressed(unsigned int key, int modifiers) override;

    inline const std::vector<std::string>& get_data_names() const
    {
        return data_names_;
    }
    std::shared_ptr<igl::opengl::ViewerDataExt>
    get_data(const std::string& data) const;

    void launch();
    void load_scene();
    void redraw_scene();
    bool pre_draw_loop();

    inline bool load(std::string scene_filename) override
    {
        if (m_state.load_scene(scene_filename)) {
            load_scene();
            return true;
        }
        return false;
    }
    inline bool save(std::string scene_filename) override
    {
        m_state.save_simulation(scene_filename);
        return true;
    }

    inline void reload()
    {
        m_reloading_scene = true;
        m_state.reload_scene();
        load_scene();
        m_reloading_scene = false;
    }

    inline void simulation_step()
    {
        // PROFILE_MAIN_POINT("simulation_step")
        // PROFILE_START()
        m_state.simulation_step();
        m_state.save_simulation_step();
        // PROFILE_END()
        // LOG_PROFILER(m_state.scene_file);
        redraw_scene();
    }

    inline void solve_collisions()
    {
        m_state.solve_collision();
        m_state.save_simulation_step();
        redraw_scene();
    }
    inline void step_solve_collisions()
    {
        m_state.collision_resolution_step();
        m_state.save_simulation_step();
        redraw_scene();
    }

    igl::opengl::glfw::Viewer m_viewer;
    SimState m_state;
    PlayerState m_player_state;

    bool m_has_scene;
    bool m_bkp_had_collision;
    bool m_bkp_has_intersections;
    int m_log_level;        ///< @brief setup log
    double m_interval_time; ///< @brief time within the interval
    bool m_show_vertex_data;

protected:
    void draw_io();
    void draw_simulation_player();
    void draw_settings();
    void draw_collision_menu();
    void draw_legends();

private:
    std::shared_ptr<igl::opengl::MeshData> mesh_data;
    std::shared_ptr<igl::opengl::VectorFieldData> velocity_data;

    std::map<std::string, std::shared_ptr<igl::opengl::ViewerDataExt>> datas_;
    std::vector<std::string> data_names_;

    bool m_reloading_scene;
};

} // namespace ccd
