#ifndef CCD_VIEWER_H
#define CCD_VIEWER_H

#include <array>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <opt/qp_solver.hpp>

#include "state.hpp"

namespace ccd {

enum ViewerEditMode {
    select,
    translate,
    add_node,
    none,
};

static const std::array<ViewerEditMode, 3> ViewerEditModeAll
    = { { select, translate, add_node } };
static const std::array<std::string, 3> ViewerEditModeNames
    = { { "Select", "Translate", "Add Node" } };

class ViewerMenu : public igl::opengl::glfw::imgui::ImGuiMenu {
private:
    typedef igl::opengl::glfw::imgui::ImGuiMenu Super;

public:
    ViewerMenu(const std::string scene_file = "");

    virtual void init(igl::opengl::glfw::Viewer* _viewer) override;
    virtual bool mouse_down(int button, int modifier) override;
    virtual bool key_pressed(unsigned int key, int modifiers) override;
    virtual void draw_menu() override;

    bool load_scene();
    bool load_scene(const std::string filename);
    bool save_scene();
    bool save_scene(const std::string filename);
    void load_state();

    // DRAW actions
    void create_edges();
    void redraw_at_time();
    void redraw_edges(const Eigen::MatrixXd& vertices);
    void recolor_edges();

    void create_displacements();
    void redraw_displacements();
    void recolor_displacements();

    void create_opt_displacements();
    void redraw_opt_displacements();
    void recolor_opt_displacements();

    void create_grad_volume();
    void redraw_grad_volume(const bool use_opt_volume);
    void recolor_grad_volume();

    void redraw_at_opt_time();
    void resize_canvas();

    void redraw_scene();

    // CRUD actions
    void connect_selected_vertices();
    void clicked__select(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void clicked__translate(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void clicked__add_node(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void undo();
    void clicked_on_canvas(
        const int button, const int modifier, const Eigen::RowVector2d& coord);

    // CCD actions
    void compute_collisions();

    // OPT actions
    void optimize_displacements();
    void load_optimization();
    void save_optimization();

    // menu windows
    void draw_io();
    void draw_edit_modes();
    void draw_ccd_steps();
    void draw_legends();
    void draw_optimization();
    void draw_optimization_results();

    State state;
    std::vector<State> state_history;

    ViewerEditMode edit_mode;

    // passive colors
    Eigen::RowVector3d color_canvas;
    Eigen::RowVector3d color_edge;
    Eigen::RowVector3d color_displ;
    Eigen::RowVector3d color_grad;
    Eigen::RowVector3d color_opt_displ;

    // active colors
    Eigen::RowVector3d color_sl;

    unsigned long canvas_data_id;
    unsigned long edges_data_id;
    unsigned long displ_data_id;
    unsigned long opt_displ_data_id;
    unsigned long gradient_data_id;
    unsigned long volume_data_id;

    Eigen::MatrixXd canvas_nodes;
    Eigen::MatrixXi canvas_faces;

    std::string scene_file;
    std::string default_scene = R"(
                                {
                                "vertices":[[-1.0, 0.0],[1.0,0.0],[0.0,0.5],[0.0,1.5]],
                                "edges":[[0,1],[2,3]],
                                "displacements":[[0.0,0.5],[0.0,0.5],[0.0,-1.0],[0.0,-1.0]]
                                })";
    std::string last_action_message;
    bool last_action_success;
};
} // namespace ccd

#endif
