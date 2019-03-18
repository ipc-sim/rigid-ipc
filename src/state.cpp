#include <algorithm> // std::sort
#include <iostream>

#include "state.hpp"

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <io/opt_results.hpp>
#include <io/read_scene.hpp>
#include <io/write_scene.hpp>

#include <opt/displacement_opt.hpp>

namespace ccd {
State::State()
    : canvas_width(10)
    , canvas_height(10)
    , current_ev_impact(-1)
    , current_edge(-1)
    , min_edge_width(0.0)
{
}

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);

    // fit scene to canvas
    Eigen::MatrixX2d all_vertices(vertices.rows() * 2, 2);
    all_vertices << vertices, vertices + displacements;

    Eigen::Vector2d v_min = all_vertices.colwise().minCoeff();
    Eigen::Vector2d v_max = all_vertices.colwise().maxCoeff();
    Eigen::RowVector2d center = all_vertices.colwise().mean();

    Eigen::Vector2d bbox = v_max - v_min;
    if (bbox[0] > canvas_width || bbox[1] > canvas_height) {
        double scale = std::min(
            canvas_width * 0.5 / bbox[0], canvas_height * 0.5 / bbox[1]);

        vertices = (vertices.rowwise() - center) * scale;
        displacements = displacements * scale;
    }

    reset_scene();
}

void State::reset_scene()
{
    reset_impacts();

    current_edge = -1;
    current_ev_impact = -1;
    time = 0.0;
    opt_time = 0.0;
    selected_displacements.clear();
    selected_points.clear();

    opt_results.x.resizeLike(displacements);
    opt_results.x.setZero();
    opt_results.minf = 2e19;
    opt_results.success = false;
}

void State::save_scene(std::string filename)
{
    io::write_scene(filename, vertices, edges, displacements);
}
// -----------------------------------------------------------------------------
// CRUD Scene
// -----------------------------------------------------------------------------
void State::add_vertex(const Eigen::RowVector2d& position)
{
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, kDIM);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, kDIM);
    displacements.row(lastid) << 0.0, -0.1;

    opt_results.x.conservativeResize(lastid + 1, kDIM);
    opt_results.x.setZero();

    reset_impacts();
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);

    // Add a new rows to the volume vector
    edge_impact_map.conservativeResize(lastid + new_edges.rows());
    volumes.conservativeResize(lastid + new_edges.rows());

    reset_impacts();
}

void State::set_vertex_position(
    const int vertex_idx, const Eigen::RowVector2d& position)
{
    vertices.row(vertex_idx) = position;
    reset_impacts();
}

void State::move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta)
{
    vertices.row(vertex_idx) += delta;
    reset_impacts();
}

void State::move_displacement(
    const int vertex_idx, const Eigen::RowVector2d& delta)
{
    displacements.row(vertex_idx) += delta;
    reset_impacts();
}

Eigen::MatrixX2d State::get_vertex_at_time()
{
    return vertices + displacements * double(time);
}

Eigen::MatrixX2d State::get_volume_grad()
{
    if (current_edge < 0 || volume_grad.cols() == 0) {
        return Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    }
    Eigen::MatrixXd grad = volume_grad.col(current_edge);
    grad.resize(grad.rows() / kDIM, kDIM);

    return grad;
}

const EdgeEdgeImpact& State::get_edge_impact(const int edge_id)
{
    size_t ee_impact_id = size_t(edge_impact_map[edge_id]);
    assert(ee_impact_id >= 0);
    return ee_impacts[ee_impact_id];
}

// -----------------------------------------------------------------------------
// CCD
// -----------------------------------------------------------------------------
void State::reset_impacts()
{

    volumes.resize(edges.rows());
    volumes.setZero();

    edge_impact_map.resize(edges.rows());
    edge_impact_map.setConstant(-1);

    volume_grad.resize(vertices.size(), edges.rows());
    volume_grad.setZero();

    ev_impacts.clear();
    ee_impacts.clear();
}

void State::detect_edge_vertex_collisions()
{
    // Get impacts between vertex and edge
    ccd::detect_edge_vertex_collisions(
        vertices, displacements, edges, ev_impacts, detection_method);

    // Sort impacts by time for convient visualization
    std::sort(ev_impacts.begin(), ev_impacts.end(),
        compare_impacts_by_time<EdgeVertexImpact>);

    // Transform to impacts between two edges
    convert_edge_vertex_to_edge_edge_impacts(
        this->edges, this->ev_impacts, this->ee_impacts);

    // Assign first impact to each edge; we store one impact for each edge on
    // edges_impact and the impacts in ee_impacts
    this->num_pruned_impacts
        = ccd::prune_impacts(this->ee_impacts, this->edge_impact_map);
}

void State::compute_collision_volumes()
{
    assert(this->volume_grad.cols() == this->edges.rows());
    assert(this->volume_grad.rows() == this->vertices.size());
    ccd::compute_volumes_fixed_toi(vertices, displacements, edges, ee_impacts,
        edge_impact_map, volume_epsilon, volumes);

    ccd::autodiff::compute_volumes_gradient(vertices, displacements, edges,
        ee_impacts, edge_impact_map, volume_epsilon, volume_grad);
}

void State::run_full_pipeline()
{
    this->detect_edge_vertex_collisions();
    this->compute_collision_volumes();
}

// -----------------------------------------------------------------------------
// OPT
// -----------------------------------------------------------------------------

void State::optimize_displacements()
{
    // reset results
    opt_results.minf = 1E10;
    opt_results.x = Eigen::MatrixXd();
    opt_results.success = false;
    opt_results.finished = false;
    u_history.clear();
    f_history.clear();
    g_history.clear();

    opt_results.x.resizeLike(displacements);
    if (!this->reuse_opt_displacements) {
        opt_results.x.setZero();
    }

    ccd::opt::setup_displacement_optimization_problem(vertices, displacements,
        edges, volume_epsilon, detection_method, opt_problem);

    Eigen::MatrixX2d U0;
    if (solver_settings.method == opt::LINEARIZED_CONSTRAINTS) {
        U0 = displacements;
    } else {
        U0 = opt_results.x;
    }

    opt_results = ccd::opt::displacement_optimization(
        opt_problem, U0, u_history, f_history, g_history, solver_settings);

    opt_results.method = solver_settings.method;
    opt_results.finished = true;

    // ui elements
    opt_time = 1.0;
    opt_iteration = -1;
}

Eigen::MatrixX2d State::get_opt_vertex_at_time(const int iteration)
{
    return vertices + get_opt_displacements(iteration) * double(opt_time);
}

Eigen::MatrixX2d State::get_opt_displacements(const int iteration)
{
    auto displ = opt_results.x;
    if (iteration >= 0 && u_history.size() > 0) {
        opt_iteration %= u_history.size();
        displ = u_history[size_t(opt_iteration)];
    }
    return displ;
}

double State::get_opt_functional(const int iteration)
{
    auto fun = opt_results.minf;
    if (iteration >= 0 && f_history.size() > 0) {
        opt_iteration %= f_history.size();
        fun = f_history[size_t(opt_iteration)];
    }
    return fun;
}

void State::save_optimization(std::string filename)
{
    io::write_opt_results(
        filename, opt_results, u_history, f_history, g_history);
}

void State::load_optimization(std::string filename)
{
    io::read_opt_results(
        filename, opt_results, u_history, f_history, g_history);

    opt_iteration = -1;
    opt_time = 1.0;
}

} // namespace ccd
