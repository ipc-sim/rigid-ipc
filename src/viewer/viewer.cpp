#include "viewer.hpp"
#include <sstream>

#include <igl/bounding_box.h>
#include <igl/opengl/MeshGL.h>

#include <ccd/not_implemented_error.hpp>

#include <io/read_scene.hpp>
#include <viewer/edges_to_rectangles.hpp>

namespace ccd {

ViewerMenu::ViewerMenu(std::string scene_file)
    : edit_mode(ViewerEditMode::none)
    , color_canvas(0.3, 0.3, 0.5)          // #4c4c80
    , color_edge(1.0, 0.0, 0.0)            // #ff0000
    , color_displ(0.0, 1.0, 0.0)           // #00ff00
    , color_grad(1.0, 0.5, 0.0)            // #ff8000
    , color_opt_displ(0.1875, 0.5, 0.1875) // #308030
    , color_collision(0.0, 1.0, 1.0)        // #00ffff
    , color_sl(1.0, 1.0, 0.0)              // #ffff00
    , scene_file(scene_file)
    , last_action_message("")
    , last_action_success(true)

{
    // clang-format off
    canvas_nodes = (Eigen::MatrixXd(4, 3) <<
        -.5, -.5, 0.0,
        -.5, 0.5, 0.0,
        0.5, -.5, 0.0,
        0.5, 0.5, 0.0).finished();

    canvas_faces = (Eigen::MatrixXi(2, 3) <<
        0, 2, 1,
        2, 3, 1).finished();
    // clang-format on
}

void ViewerMenu::init(igl::opengl::glfw::Viewer* _viewer)
{
    Super::init(_viewer);

    // first added is drawn in front
    gradient_data_id = viewer->data_list.size() - 1;

    viewer->append_mesh();
    edges_data_id = viewer->data_list.size() - 1;
    viewer->data().show_vertid = true;

    viewer->append_mesh();
    opt_displ_data_id = viewer->data_list.size() - 1;

    viewer->append_mesh();
    displ_data_id = viewer->data_list.size() - 1;

    viewer->append_mesh();
    volume_data_id = viewer->data_list.size() - 1;

    load_scene(scene_file);

    // CANVAS
    {
        viewer->append_mesh();
        viewer->data().set_mesh(canvas_nodes, canvas_faces);
        viewer->data().show_lines = false;
        viewer->data().set_colors(color_canvas);
        canvas_data_id = viewer->data_list.size() - 1;
        resize_canvas();
    }
}

bool ViewerMenu::save_scene()
{
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
        return false;
    state.save_scene(fname);
    return true;
}

bool ViewerMenu::load_scene()
{
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
        return false;
    return load_scene(fname);
}

bool ViewerMenu::load_scene(const std::string filename)
{

    if (filename.empty()) {
        io::read_scene_from_str(
            default_scene, state.vertices, state.edges, state.displacements);
        state.reset_scene();

    } else {
        state.load_scene(filename);
    }

    state_history.clear();
    state_history.push_back(state);

    load_state();
    return true;
}

void ViewerMenu::load_state()
{
    create_edges();
    create_displacements();
    create_opt_displacements();
    recolor_opt_collisions();
    create_grad_volume();
    redraw_at_time();

}

void ViewerMenu::resize_canvas()
{
    Eigen::MatrixXd V = canvas_nodes;
    V.col(0) *= state.canvas_width;
    V.col(1) *= state.canvas_height;
    viewer->data_list[canvas_data_id].set_vertices(V);
}

// ----------------------------------------------------------------------------------------------------------------------------
// CCD USER ACTIONS
// ----------------------------------------------------------------------------------------------------------------------------
void ViewerMenu::compute_collisions()
{
    try {
        state.run_ccd_pipeline();
        redraw_grad_volume(/*use_opt_grad=*/false);
        state_history.push_back(state);

    } catch (NotImplementedError e) {
        last_action_message = e.what();
        last_action_success = false;
    }
}


// ----------------------------------------------------------------------------------------------------------------------------
// OPT USER ACTIONS
// ----------------------------------------------------------------------------------------------------------------------------
void ViewerMenu::optimize_displacements()
{
    try {
        state.optimize_displacements();
        last_action_success = state.opt_results.success;
        std::ostringstream message;
        message.precision(3);
        message << "Optimization " << (last_action_success ? "" : "un")
                << "successful" << std::endl;
        this->last_action_message = message.str();
        redraw_opt_displacements();
        redraw_at_opt_time();
    } catch (NotImplementedError e) {
        last_action_message = e.what();
        last_action_success = false;
    } catch (std::runtime_error e) {
        last_action_message = e.what();
        last_action_success = false;
    }
}

void ViewerMenu::load_optimization()
{
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0) {
        return;
    }
    state.load_optimization(fname);
    redraw_opt_displacements();
    redraw_at_opt_time();
}

void ViewerMenu::save_optimization()
{
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0) {
        return;
    }
    return state.save_optimization(fname);
}

void ViewerMenu::post_process_optimization(){
    state.post_process_optimization();
}

} // namespace ccd
