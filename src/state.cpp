#include <algorithm> // std::sort
#include <iostream>

#include <igl/slice.h>

#include "state.hpp"

#include <ccd/collision_penalty_diff.hpp>
#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <io/opt_results.hpp>
#include <io/read_scene.hpp>
#include <io/write_scene.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/barrier_constraint.hpp>
#include <opt/displacement_opt.hpp>
#include <opt/ncp_solver.hpp>
#include <opt/volume_constraint.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

State::State()
    : convert_to_rigid_bodies(false)
    , output_dir(DATA_OUTPUT_DIR)
    , opt_method(OptimizationMethod::BARRIER_NEWTON)
    , constraint_function(ConstraintType::BARRIER)
    , log_level(spdlog::level::info)
    , canvas_width(10)
    , canvas_height(10)
    , current_volume(-1)
    , grad_scaling(1.0f)
    , use_opt_gradient(false)
    , current_opt_time(0.0)
    , current_opt_iteration(-1)
{
    barrier_newton_solver.barrier_constraint = &barrier_constraint;
    opt_problem.intermediate_callback
        = [&](const Eigen::VectorXd& x, const Eigen::MatrixX2d& Uk) {
              return record_optimization_step(x, Uk);
          };

    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
}

// -----------------------------------------------------------------------------
// CRUD Scene
// -----------------------------------------------------------------------------

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);
    if (convert_to_rigid_bodies) {
        convert_connected_components_to_rigid_bodies();
    }
    fit_scene_to_canvas();
    reset_scene();
}

void State::load_scene(const Eigen::MatrixX2d& v, const Eigen::MatrixX2i& e,
    const Eigen::MatrixX2d& d)
{
    this->vertices = v;
    this->edges = e;
    this->displacements = d;

    reset_scene();
}

void State::fit_scene_to_canvas()
{

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
}

void State::reset_scene()
{
    reset_scene_data();

    current_volume = -1;
    current_time = 0.0;
    current_opt_time = 0.0;
    selected_displacements.clear();
    selected_points.clear();

    opt_results.x.resizeLike(displacements);
    opt_results.x.setZero();
    opt_results.minf = 2e19;
    opt_results.success = false;
    opt_results.finished = false;

    opt_problem.fixed_dof.resize(displacements.size(), 1);
    opt_problem.fixed_dof.setConstant(false);
}

void State::save_scene(std::string filename)
{
    io::write_scene(filename, vertices, edges, displacements);
}

void State::add_vertex(const Eigen::RowVector2d& position)
{
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, kDIM);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, kDIM);
    displacements.row(lastid) << 0.0, -0.1;

    opt_results.x.conservativeResize(lastid + 1, kDIM);
    opt_results.x.setZero();

    reset_scene_data();

    // Resize the fixed_dof
    expand_fixed_dof();
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);

    reset_scene_data();
}

void State::expand_fixed_dof()
{
    long num_old_vertices = opt_problem.fixed_dof.size() / 2;
    long num_new_vertices = displacements.rows() - num_old_vertices;

    opt_problem.fixed_dof.resize(num_old_vertices, 2);
    opt_problem.fixed_dof.conservativeResize(displacements.rows(), 2);
    opt_problem.fixed_dof.bottomRows(num_new_vertices).setZero();
    opt_problem.fixed_dof.resize(opt_problem.fixed_dof.size(), 1);

    assert(opt_problem.fixed_dof.size() == displacements.size());
}

// Remove vertices that are not an end-point to any edge.
bool State::remove_free_vertices()
{
    // Create a set of non-free vertices to keep
    std::set<int> nonfree_ids_set;
    for (int i = 0; i < edges.rows(); i++) {
        nonfree_ids_set.insert(edges(i, 0));
        nonfree_ids_set.insert(edges(i, 1));
    }
    std::vector<int> nonfree_ids
        = std::vector<int>(nonfree_ids_set.begin(), nonfree_ids_set.end());
    if (nonfree_ids.size() == vertices.rows()) {
        return false; // Do nothing
    }

    // Convert the set of non-free vertices to an Eigen vector
    Eigen::VectorXi R
        = Eigen::VectorXi::Map(nonfree_ids.data(), long(nonfree_ids.size()));

    // Create a range for the columns to keep
    Eigen::VectorXi C = Eigen::VectorXi::LinSpaced(
        displacements.cols(), 0, displacements.cols() - 1);

    // Create a temporary matrix of the vertices and displacements to keep
    Eigen::MatrixX2d new_vertices, new_displacements;
    Eigen::MatrixXb new_fixed_dof;
    igl::slice(vertices, R, C, new_vertices); // Copy over the vertices to keep
    igl::slice(displacements, R, C, new_displacements); // Copy displacements
    opt_problem.fixed_dof.resize(opt_problem.fixed_dof.size() / 2, 2);
    igl::slice(opt_problem.fixed_dof, R, C, new_fixed_dof);
    vertices = new_vertices;
    displacements = new_displacements;

    // Update the edges to have the new vertex indices
    // Map from old to new indices
    std::unordered_map<int, int> old_ids_to_new_ids;
    int new_id = 0;
    for (const auto& old_id : nonfree_ids) {
        // The new indices are the order of indices in nonfree_vertex_ids
        old_ids_to_new_ids.insert(std::make_pair(old_id, new_id++));
    }
    for (int i = 0; i < edges.rows(); i++) {
        edges(i, 0) = old_ids_to_new_ids.at(edges(i, 0));
        edges(i, 1) = old_ids_to_new_ids.at(edges(i, 1));
    }

    reset_scene();
    new_fixed_dof.resize(new_fixed_dof.size(), 1);
    opt_problem.fixed_dof = new_fixed_dof;

    return true;
}

// Duplicate selected vertices and edges that have both end-points selected.
void State::duplicate_selected_vertices(Eigen::Vector2d delta_center_of_mass)
{
    std::vector<long> selected_edges;
    for (int i = 0; i < edges.rows(); i++) {
        if (std::find(
                selected_points.begin(), selected_points.end(), edges(i, 0))
                != selected_points.end()
            && std::find(
                   selected_points.begin(), selected_points.end(), edges(i, 1))
                != selected_points.end()) {
            selected_edges.push_back(i);
        }
    }

    long num_old_vertices = vertices.rows();
    long num_new_vertices = long(selected_points.size());

    long num_old_edges = edges.rows();
    long num_new_edges = long(selected_edges.size());

    edges.conservativeResize(num_old_edges + num_new_edges, 2);
    // Map from old to new indices
    std::unordered_map<int, int> old_ids_to_new_ids;
    long new_id = num_old_vertices;
    for (const int& selected_point : selected_points) {
        old_ids_to_new_ids.insert(std::make_pair(selected_point, new_id++));
    }
    long new_edge_id = num_old_edges;
    for (const long& selected_edge : selected_edges) {
        edges(new_edge_id, 0) = old_ids_to_new_ids.at(edges(selected_edge, 0));
        edges(new_edge_id, 1) = old_ids_to_new_ids.at(edges(selected_edge, 1));
        new_edge_id++;
    }

    vertices.conservativeResize(num_old_vertices + num_new_vertices, 2);
    displacements.conservativeResize(num_old_vertices + num_new_vertices, 2);
    Eigen::VectorXi R = Eigen::VectorXi::Map(
        selected_points.data(), long(selected_points.size()));
    Eigen::VectorXi C = Eigen::VectorXi::LinSpaced(
        displacements.cols(), 0, displacements.cols() - 1);

    Eigen::MatrixX2d duplicate_vertices, duplicate_displacements;
    igl::slice(vertices, R, C, duplicate_vertices);
    duplicate_vertices.col(0).array() += delta_center_of_mass.x();
    duplicate_vertices.col(1).array() += delta_center_of_mass.y();
    vertices.bottomRows(num_new_vertices) = duplicate_vertices;
    igl::slice(displacements, R, C, duplicate_displacements);
    displacements.bottomRows(num_new_vertices) = duplicate_displacements;

    opt_results.x.resize(num_old_vertices + num_new_vertices, 2);
    opt_results.x.setZero();

    reset_scene_data();

    // Resize the fixed_dof
    expand_fixed_dof();
}

void State::set_vertex_position(
    const int vertex_idx, const Eigen::RowVector2d& position)
{
    vertices.row(vertex_idx) = position;
    reset_scene_data();
}

void State::move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta)
{
    vertices.row(vertex_idx) += delta;
    reset_scene_data();
}

void State::move_displacement(
    const int vertex_idx, const Eigen::RowVector2d& delta)
{
    displacements.row(vertex_idx) += delta;
    reset_scene_data();
}

// -----------------------------------------------------------------------------
// CCD
// -----------------------------------------------------------------------------
void State::reset_scene_data()
{
    volumes.resize(0);
    volumes.setZero();

    volume_grad.resize(0, 0);
    volume_grad.setZero();

    volume_constraint.initialize(vertices, edges, displacements);
    barrier_constraint.initialize(vertices, edges, displacements);

    u_history.clear();
    f_history.clear();
    gsum_history.clear();
    jac_g_history.clear();
    collision_history.clear();
}

void State::run_ccd_pipeline()
{
    // NOTE: always uses volume constraint
    volume_constraint.initialize(vertices, edges, displacements);
    volume_constraint.compute_constraints(displacements, volumes);
    volume_constraint.compute_constraints_jacobian(displacements, volume_grad);

    current_volume = 0;

    const Eigen::IOFormat CommaFmt(
        Eigen::StreamPrecision, Eigen::DontAlignCols, ",", ",", "", "", "", "");

    spdlog::debug("grad_vol\n{}", ccd::log::fmt_eigen(volume_grad, 4));
}

// -----------------------------------------------------------------------------
// OPT
// -----------------------------------------------------------------------------
opt::CollisionConstraint& State::getCollisionConstraint()
{
    if (constraint_function == ConstraintType::VOLUME) {
        return volume_constraint;
    } else {
        return barrier_constraint;
    }
}

opt::OptimizationSolver& State::getOptimizationSolver()
{
    switch (opt_method) {
    case OptimizationMethod::NCP:
        return ncp_solver;
    case OptimizationMethod::IPOPT:
        return ipopt_solver;
    case OptimizationMethod::NLOPT:
        return nlopt_solver;
    case OptimizationMethod::LINEARIZED_CONSTRAINTS:
        return qp_solver;
    case OptimizationMethod::BARRIER_NEWTON:
        return barrier_newton_solver;
    }
}

void State::reset_optimization_problem()
{
    volume_constraint.initialize(vertices, edges, displacements);

    auto& constraint = getCollisionConstraint();
    constraint.initialize(vertices, edges, displacements);
    opt_problem.initialize(vertices, edges, displacements, constraint);
}

void State::reset_results()
{
    opt_results.success = false;
    opt_results.finished = false;
    bool dirty = opt_results.x.size() == 0
        || opt_results.x.rows() != displacements.rows();

    if (dirty || !reuse_opt_displacements) {
        u_history.clear();
        f_history.clear();
        gsum_history.clear();
        jac_g_history.clear();
        opt_results.minf = 1E10;
        opt_results.x.resizeLike(displacements);
    }
}

bool State::record_optimization_step(
    const Eigen::VectorXd& x, const Eigen::MatrixX2d& Uk)
{
    u_history.push_back(Uk);
    // check if there is any reason to stop optimization
    return true;
}

void State::post_process_optimization()
{
    f_history.clear();
    f_history.reserve(u_history.size());
    g_history.clear();
    g_history.reserve(u_history.size());
    gsum_history.clear();
    gsum_history.reserve(u_history.size());
    jac_g_history.clear();
    jac_g_history.reserve(u_history.size());
    collision_history.clear();

    size_t i = 0, steps = u_history.size();
    for (auto& uk : u_history) {
        Eigen::MatrixXd x = uk;
        x.resize(x.size(), 1);
        f_history.push_back(opt_problem.eval_f(x));

        Eigen::VectorXd gx;
        Eigen::MatrixXd jac_gx;
        volume_constraint.detectCollisions(uk);
        volume_constraint.compute_active_constraints(uk, gx, jac_gx);
        g_history.push_back(gx);
        gsum_history.push_back(gx.sum());
        jac_g_history.push_back(jac_gx);

        std::vector<long> active_edges;
        active_edges.reserve(volume_constraint.ee_impacts.size());
        for (auto& ee : volume_constraint.ee_impacts) {
            active_edges.push_back(ee.impacted_edge_index);
        }
        collision_history.push_back(active_edges);

        spdlog::debug("post-process {}/{} gx_sum={}", ++i, steps, gx.sum());
    }
}
void State::optimize_displacements(const std::string /*filename*/)
{
    reset_optimization_problem();
    reset_results();

    // 3. set initial value
    Eigen::MatrixXd U0(displacements.rows(), 2);
    if (opt_method == OptimizationMethod::LINEARIZED_CONSTRAINTS) {
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

    if (opt_problem.validate_problem()) {
        Eigen::MatrixXd x0 = U0;
        x0.resize(U0.size(), 1);
        opt_problem.x0 = x0;

        auto& solver = getOptimizationSolver();
#ifdef PROFILE_FUNCTIONS
        reset_profiler();
        igl::Timer timer;
        timer.start();
#endif
        opt_results = solver.solve(opt_problem);
        opt_results.x.resize(U0.rows(), 2);

#ifdef PROFILE_FUNCTIONS
        timer.stop();
        print_profile(timer.getElapsedTime());
#endif
    }

    // log_optimization_steps(filename, it_x, it_lambda, it_gamma);
    opt_results.finished = true;

    // update ui elements
    current_opt_time = 1.0;
    current_opt_iteration = -1;
}

void State::save_optimization(std::string filename)
{
    io::write_opt_results(
        filename, opt_results, u_history, f_history, gsum_history);
}

void State::load_optimization(std::string filename)
{
    io::read_opt_results(
        filename, opt_results, u_history, f_history, gsum_history);

    current_opt_iteration = -1;
    current_opt_time = 1.0;
}

// -----------------------------------------------------------------------------
// UI Getters
// -----------------------------------------------------------------------------
Eigen::MatrixX2d State::get_vertex_at_time()
{
    return vertices + displacements * double(current_time);
}

Eigen::MatrixX2d State::get_volume_grad()
{
    if (current_volume < 0 || volume_grad.rows() == 0) {
        return Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    }
    current_volume %= volume_grad.rows();
    Eigen::MatrixXd grad = volume_grad.row(current_volume).transpose();
    grad.resize(grad.rows() / kDIM, kDIM);

    return grad;
}

Eigen::MatrixX2d State::get_opt_vertex_at_time()
{
    return vertices + get_opt_displacements() * double(current_opt_time);
}

Eigen::MatrixX2d State::get_opt_displacements()
{
    auto displ = opt_results.x;
    current_opt_iteration = std::max(-1, current_opt_iteration);
    current_opt_iteration
        = std::min(int(u_history.size()), current_opt_iteration);

    if (current_opt_iteration >= 0 && u_history.size() > 0) {
        displ = u_history[size_t(current_opt_iteration)];
    }
    return displ;
}

Eigen::MatrixX2d State::get_opt_volume_grad()
{
    current_volume = std::max(-1, current_volume);

    Eigen::MatrixX2d empty_grad = Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    if (current_opt_iteration < 0 || jac_g_history.size() == 0
        || current_volume < 0) {
        return empty_grad;
    }

    auto& grad_volume = jac_g_history[size_t(current_opt_iteration)];
    if (grad_volume.size() == 0) {
        return empty_grad;
    }

    current_volume = std::min(int(grad_volume.rows()) - 1, current_volume);

    Eigen::MatrixXd grad = grad_volume.row(current_volume).transpose();
    grad.resize(grad.rows() / kDIM, kDIM);

    return grad;
}

Eigen::VectorXd State::get_opt_volume()
{
    if (current_opt_iteration < 0 || gsum_history.size() == 0) {
        return Eigen::VectorXd::Zero(edges.rows());
    } else {
        return g_history[size_t(current_opt_iteration)];
    }
}

double State::get_opt_functional()
{
    if (current_opt_iteration == -1) {
        return opt_results.minf;
    }
    current_opt_iteration = std::max(-1, current_opt_iteration);
    if (f_history.size() > 0) {
        current_opt_iteration
            = std::min(current_opt_iteration, int(f_history.size()));
        return f_history[size_t(current_opt_iteration)];
    }
    return -1;
}

std::vector<std::vector<int>> State::create_adjacency_list()
{
    std::vector<std::vector<int>> adjacency_list(vertices.rows());
    for (int i = 0; i < edges.rows(); i++) {
        adjacency_list[edges(i, 0)].push_back(edges(i, 1));
        adjacency_list[edges(i, 1)].push_back(edges(i, 0));
    }
    return adjacency_list;
}

void State::find_connected_vertices(const long& vertex_id,
    const std::vector<std::vector<int>>& adjacency_list,
    std::unordered_set<int>& connected_vertices)
{
    if (connected_vertices.find(vertex_id) != connected_vertices.end()) {
        return;
    }
    connected_vertices.insert(vertex_id);
    for (const int& adjacent_vertex_id : adjacency_list[vertex_id]) {
        find_connected_vertices(
            adjacent_vertex_id, adjacency_list, connected_vertices);
    }
}

void State::convert_connected_components_to_rigid_bodies()
{
    auto adjacency_list = create_adjacency_list();
    rigid_bodies.clear();
    Eigen::VectorXb is_part_of_body = Eigen::VectorXb::Zero(vertices.rows());

    std::unordered_set<int> connected_vertices;
    Eigen::MatrixXi body_edges;
    Eigen::MatrixX2d body_vertices;

    for (long i = 0; i < vertices.rows(); i++) {
        if (!is_part_of_body(i)) {
            find_connected_vertices(i, adjacency_list, connected_vertices);

            std::vector<int> ordered_connected_vertices(
                connected_vertices.begin(), connected_vertices.end());
            std::unordered_map<int, int> old_ids_to_new_ids;
            int new_id = 0;
            for (const auto& old_id : ordered_connected_vertices) {
                // The new indices are the order of indices
                old_ids_to_new_ids.insert(std::make_pair(old_id, new_id++));
            }

            std::vector<int> column0, column1;
            for (long j = 0; j < edges.rows(); j++) {
                if (connected_vertices.find(edges(j, 0))
                    != connected_vertices.end()) {
                    column0.push_back(old_ids_to_new_ids.at(edges(j, 0)));
                    column1.push_back(old_ids_to_new_ids.at(edges(j, 1)));
                }
            }
            body_edges = Eigen::MatrixXi(long(column0.size()), 2);
            body_edges.col(0)
                = Eigen::VectorXi::Map(column0.data(), long(column0.size()));
            body_edges.col(1)
                = Eigen::VectorXi::Map(column1.data(), long(column1.size()));

            Eigen::VectorXi R
                = Eigen::VectorXi::Map(ordered_connected_vertices.data(),
                    long(ordered_connected_vertices.size()));
            igl::slice(vertices, R, Eigen::VectorXi::LinSpaced(2, 0, 1),
                body_vertices);

            rigid_bodies.push_back(RigidBody(
                body_vertices, body_edges, Eigen::Vector3d(0, 0, 3.14 / 4)));

            for (const auto& vertex_id : ordered_connected_vertices) {
                is_part_of_body(vertex_id) = true;
            }
            connected_vertices.clear();
        }
    }

    vertices = Eigen::MatrixX2d();
    displacements = Eigen::MatrixX2d();
    edges = Eigen::MatrixX2i();
    for (const auto& rigid_body : rigid_bodies) {
        long prev_vertex_count = vertices.rows();
        vertices.conservativeResize(
            vertices.rows() + rigid_body.vertices.rows(), 2);
        displacements.conservativeResize(
            displacements.rows() + rigid_body.vertices.rows(), 2);
        edges.conservativeResize(edges.rows() + rigid_body.edges.rows(), 2);
        vertices.bottomRows(rigid_body.vertices.rows()) = rigid_body.vertices;
        Eigen::MatrixX2d body_displacements;
        rigid_body.compute_particle_displacments(body_displacements);
        displacements.bottomRows(rigid_body.vertices.rows())
            = body_displacements;
        edges.bottomRows(rigid_body.edges.rows()) = rigid_body.edges;
        edges.bottomRows(rigid_body.edges.rows()).array() += prev_vertex_count;
    }

    reset_scene();
}

std::vector<long> State::get_opt_collision_edges()
{
    if (current_opt_iteration == -1) {
        return std::vector<long>();
    }
    current_opt_iteration = std::max(-1, current_opt_iteration);
    if (collision_history.size() > 0) {
        current_opt_iteration
            = std::min(current_opt_iteration, int(collision_history.size()));
        return collision_history[size_t(current_opt_iteration)];
    }
    return std::vector<long>();
}

} // namespace ccd
