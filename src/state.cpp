#include <algorithm> // std::sort
#include <iostream>

#include "state.hpp"

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <io/opt_results.hpp>
#include <io/read_scene.hpp>
#include <io/write_scene.hpp>

#include <logger.hpp>

#include <opt/displacement_opt.hpp>

namespace ccd {
State::State()
    : volume_epsilon(1E-3)
    , recompute_collision_set(false)
    , canvas_width(10)
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
    opt_results.finished = false;
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
    // 1. setup problems
    ccd::opt::setup_displacement_optimization_problem(vertices, displacements,
        edges, volume_epsilon, detection_method, recompute_collision_set,
        opt_problem);

    // 2. reset results
    opt_results.success = false;
    opt_results.finished = false;
    bool dirty = opt_results.x.size() == 0
        || opt_results.x.rows() != displacements.rows();
    if (dirty || !reuse_opt_displacements) {
        u_history.clear();
        f_history.clear();
        g_history.clear();
        opt_results.minf = 1E10;
        opt_results.x.resizeLike(displacements);
    }

    // 3. set initial value
    Eigen::MatrixXd U0(displacements.rows(), 2);
    if (solver_settings.method == opt::LINEARIZED_CONSTRAINTS) {
        U0 = displacements;
    } else if (reuse_opt_displacements) {
        U0 = opt_results.x;
    } else {
        U0.setZero();
    }
    opt_results.x = U0;

    // 4. setup callback
    std::vector<Eigen::VectorXd> it_x;
    std::vector<Eigen::VectorXd> it_lambda;
    std::vector<double> it_gamma;

    solver_settings.intermediate_cb
        = [&](const Eigen::VectorXd& x, const double obj_value,
              const Eigen::VectorXd& dual, const double gamma,
              const int /*iteration*/) {
              Eigen::MatrixXd u = x;
              u.resize(displacements.rows() / 2, 2);
              u_history.push_back(u);
              f_history.push_back(obj_value);
              auto gx = opt_problem.g(x);
              g_history.push_back(gx.sum());

              it_x.push_back(x);
              it_lambda.push_back(dual);
              it_gamma.push_back(gamma);
          };
    // 5. run optimization
    opt_results
        = ccd::opt::displacement_optimization(opt_problem, U0, solver_settings);

    // LOGGING for testing only --- -should move somewhere else
    const Eigen::IOFormat CommaFmt(Eigen::StreamPrecision,
        Eigen::DontAlignCols, ", ", ", ", "", "", "", "");
    const std::string base_dir = TEST_OUTPUT_DIR;
    std::string uuid = ccd::log::now();
    std::string file_name = fmt::format("opt_{0}.csv", uuid);
    std::ofstream o(base_dir + "/" + file_name);

    // print table header
    o << ccd::opt::OptimizationMethodNames[solver_settings.method] << ",";
    o << ccd::opt::LCPSolverNames[solver_settings.lcp_solver] << std::endl;
    o << "it";
    for (uint i=0; i < opt_problem.num_vars; i++){
        o << ",x_" << i ;
    }
    for (uint i=0; i < opt_problem.num_constraints; i++){
        o << ",lambda_" << i;
    }
    o << ",gamma, f(x)";
    for (uint i=0; i < opt_problem.num_constraints; i++){
        o << ",gx_" << i;
    }
    o << ",||jac_g(x)||"  << std::endl;

    // print initial value
    Eigen::MatrixXd d = displacements;
    d.resize(d.size(), 1);
    o << 0 << ",";
    o << d.format(CommaFmt) << ",";
    o << Eigen::VectorXd::Zero(opt_problem.num_constraints).format(CommaFmt) << ",";
    o << 0.0 << ",";
    o << opt_problem.f(d) << ",";
    o << opt_problem.g(d).format(CommaFmt) << ",";
    o << opt_problem.jac_g(d).squaredNorm() << std::endl << std::endl;

    // print steps
    for (uint i = 0; i < it_x.size(); i++) {
        o << i + 1 << ",";
        o << it_x[i].format(CommaFmt) << ",";
        o << it_lambda[i].format(CommaFmt) << ",";
        o << it_gamma[i] << ",";
        o << opt_problem.f(it_x[i]) << ",";
        o << opt_problem.g(it_x[i]).format(CommaFmt) << ",";
        o << opt_problem.jac_g(it_x[i]).squaredNorm() << std::endl;
    }

    opt_results.method = solver_settings.method;
    opt_results.finished = true;

    // update ui elements
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
