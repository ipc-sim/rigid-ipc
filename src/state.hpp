#ifndef CCD_STATE_HPP
#define CCD_STATE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <FixingCollisions/collision_detection.hpp>

namespace ccd {
class State {
public:
    ~State() = default;
    State();

    Eigen::MatrixXd vertices;
    Eigen::MatrixX2i edges;
    Eigen::MatrixXd displacements;
    Eigen::MatrixXd volume_grad;

    ImpactsPtr impacts;
    DetectionMethod detection_method;

    // ----------------------------------- SCENE CRUD
    void load_scene(const std::string filename);
    void save_scene(const std::string filename);

    void add_vertex(const Eigen::RowVector3d& vertex);
    void add_edges(const Eigen::MatrixX2i& edges);

    void set_vertex_position(const int vertex_idx, const Eigen::RowVector3d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector3d& delta);
    void move_displacement(const int vertex_idx, const Eigen::RowVector3d& delta);

    void set_vertex_displacement(const int vertex_idx, const Eigen::Vector2d& displacement);

    void detect_edge_vertex_collisions();

    // ----------------------------------- UI
    double canvas_width, canvas_height;
    float time;
    std::vector<int> selected_points, selected_displacements;
    int current_impact;


};

}
#endif
