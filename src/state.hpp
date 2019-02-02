#ifndef CCD_STATE_HPP
#define CCD_STATE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ccd {
class State {
public:
    ~State() = default;
    State();

    Eigen::MatrixXd vertices;
    Eigen::MatrixXi edges;
    Eigen::MatrixXd displacements;
    Eigen::MatrixXd volume_grad;

    // ----------------------------------- SCENE CRUD
    void load_scene(const std::string filename);

    void add_vertex(const Eigen::RowVector3d& vertex);
    void add_edges(const Eigen::MatrixXi& edges);

    void set_vertex_position(const int vertex_idx, const Eigen::RowVector3d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector3d& delta);
    void move_displacement(const int vertex_idx, const Eigen::RowVector3d& delta);

    void set_vertex_displacement(const int vertex_idx, const Eigen::Vector2d& displacement);

    // ----------------------------------- UI
    double canvas_width, canvas_height;
    float displacement_ptge;
    std::vector<int> selected_points, selected_displacements;


};

}
#endif
