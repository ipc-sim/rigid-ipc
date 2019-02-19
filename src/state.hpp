#ifndef CCD_STATE_HPP
#define CCD_STATE_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ccd/collision_detection.hpp>
#include <ccd/impact.hpp>

namespace ccd {

class State {
public:
    static const int kDIM = 2;

    State();
    virtual ~State() = default;

    Eigen::MatrixX2d vertices;
    Eigen::MatrixX2i edges;
    Eigen::MatrixX2d displacements;
    Eigen::VectorXd volumes;
    Eigen::MatrixX2d volume_grad;

    EdgeVertexImpactsPtr impacts;
    DetectionMethod detection_method;
    double epsilon = 0.0;

    // ----------------------------------- SCENE CRUD
    void load_scene(const std::string filename);
    void save_scene(const std::string filename);

    void add_vertex(const Eigen::RowVector2d& vertex);
    void add_edges(const Eigen::MatrixX2i& edges);

    void set_vertex_position(
        const int vertex_idx, const Eigen::RowVector2d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta);
    void move_displacement(
        const int vertex_idx, const Eigen::RowVector2d& delta);

    // ----------------------------------- SCENE CCD
    void detect_edge_vertex_collisions();
    void compute_collision_volumes();
    void run_full_pipeline();

    // --------------------------------------- UI
    // Background rectangle to detect clicks
    double canvas_width, canvas_height;
    // We show the scene at time=`time` between 0 and 1
    float time;
    // Current user-selection of vertex and displacement points
    std::vector<int> selected_points, selected_displacements;
    // Use for any functionallity that requires showing only one impact
    // e.g goto time_of_impact
    int current_impact;
};

}
#endif
