#pragma once

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <viewer/igl_viewer_ext.hpp>

#include "SimState.hpp"

namespace ccd {

class UISimState : public igl::opengl::glfw::imgui::ImGuiMenu {
    typedef igl::opengl::glfw::imgui::ImGuiMenu Super;

private:
    std::unique_ptr<igl::opengl::ViewerDataExt> edges_data;
    enum PlayerState { Playing = 0, Paused, TotalPlayerStatus };

public:
    UISimState();
    ~UISimState() override {}

    virtual void init(igl::opengl::glfw::Viewer* _viewer) override;
    virtual void draw_menu() override;
    //    virtual bool mouse_down(int button, int modifier) override;
    //    virtual bool key_pressed(unsigned int key, int modifiers) override;

    void launch();
    void load_scene();
    void redraw_scene();
    bool pre_draw_loop();

    igl::opengl::glfw::Viewer m_viewer;
    SimState m_state;
    PlayerState m_player_state;
    bool m_has_scene;
    bool m_bkp_had_collision;
    bool m_bkp_has_collision;

    int m_log_level; ///< brief setup log

    Eigen::RowVector3d color_mesh;

    bool load(std::string scene_filename) override
    {
        m_state.load_scene(scene_filename);
        load_scene();
        return true;
    }

    void reload(){
        m_state.reload_scene();
        load_scene();
    }

    void simulation_step()
    {
        m_state.simulation_step();
        redraw_scene();
    }

protected:
    void draw_io();
    void draw_simulation_player();

};

} // namespace ccd
