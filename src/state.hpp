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

    Eigen::MatrixX2d vertices;
    Eigen::MatrixX2i edges;
    Eigen::MatrixX2d displacements;
    Eigen::MatrixX2d volume_grad;

    ImpactsPtr impacts;
    DetectionMethod detection_method;

    // ----------------------------------- SCENE CRUD
    void load_scene(const std::string filename);
    void save_scene(const std::string filename);

    void add_vertex(const Eigen::RowVector2d& vertex);
    void add_edges(const Eigen::MatrixX2i& edges);

    void set_vertex_position(const int vertex_idx, const Eigen::RowVector2d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta);
    void move_displacement(const int vertex_idx, const Eigen::RowVector2d& delta);

        // ----------------------------------- SCENE CCD
    void detect_edge_vertex_collisions();

    // --------------------------------------- UI
    // bacground rectangle to detect clicks
    double canvas_width, canvas_height;
    // we show the scene at time=`time` between 0 and 1
    float time;
    // current user-selection of vertex and displacement points
    std::vector<int> selected_points, selected_displacements;
    // use for any functionallity that requires showing only one impact
    // e.g goto time_of_impact
    int current_impact;

};

}
#endif
