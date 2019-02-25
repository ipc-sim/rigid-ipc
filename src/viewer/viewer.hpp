#ifndef CCD_VIEWER_H
#define CCD_VIEWER_H

#include <array>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include "state.hpp"

namespace ccd {

enum ViewerEditMode {
    select,
    translate,
    add_node,
    add_chain,
};

static const std::array<ViewerEditMode, 4> ViewerEditModeAll
    = { { select, translate, add_node, add_chain } };
static const std::array<std::string, 4> ViewerEditModeNames
    = { { "Select", "Translate", "Add Node", "Add Chain" } };

static const std::array<DetectionMethod, 2> DetectionMethodAll
    = { { BRUTE_FORCE } };
static const std::array<std::string, 2> DetectionMethodNames
    = { { "BRUTE_FORCE" } };

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

    void resize_canvas();
    void clicked_on_canvas(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void recolor_vertices();
    void redraw_scene();
    void redraw_displacements();
    void redraw_at_time();
    void redraw_volumes();
    void redraw_volumes_grad();

    // CRUD actions
    void connect_selected_vertices();
    void clicked__select(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void clicked__translate(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void clicked__add_node(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void clicked__add_chain(
        const int button, const int modifier, const Eigen::RowVector2d& coord);
    void undo();

    // CCD actions
    void detect_edge_vertex_collisions();
    void compute_collision_volumes();
    void goto_ev_impact(const int ev_impact);
    void goto_ee_impact(const int ee_impact);

    // menu windows
    void draw_io();
    void draw_edit_modes();
    void draw_ui_settings();
    void draw_ccd_steps();

    // utils functions
    void update_vector_field(const unsigned long data_id,
        const Eigen::MatrixXd& x0, const Eigen::MatrixXd& delta);

    void extend_vector_field(const unsigned long data_id,
        const Eigen::MatrixXd& x0, const Eigen::MatrixXd& delta, const int last,
        const Eigen::RowVector3d& color);

    void update_graph(const unsigned long data_id, const Eigen::MatrixXd& nodes,
        const Eigen::MatrixXi& edges);

    void add_graph_vertex(const unsigned long data_id,
        const Eigen::MatrixXd& vertices, const Eigen::RowVector2d& vertex,
        const Eigen::RowVector3d& color);

    void add_graph_edges(const unsigned long data_id,
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& new_edges,
        const Eigen::RowVector3d& color);

    void color_points(
        const unsigned long data_id, const Eigen::RowVector3d& color);

    void highlight_points(const unsigned long data_id,
        const std::vector<int>& nodes, const Eigen::RowVector3d& color_hl);

    // fixed on libigl viewer functions
    void viewer_set_edges(const unsigned long data_id, const Eigen::MatrixXd& P,
        const Eigen::MatrixXi& E, const Eigen::MatrixXd& C);

    void viewer_set_vertices(
        const unsigned long data_id, const Eigen::MatrixXd& V);

    State state;
    std::vector<State> state_history;

    ViewerEditMode edit_mode;

    Eigen::RowVector3d color_vtx, color_edge, color_displ, color_grad,
        color_canvas, color_sl;
    unsigned long canvas_data_id;
    unsigned long surface_data_id;
    unsigned long displ_data_id;
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
}

#endif
