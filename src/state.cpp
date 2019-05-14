#include <algorithm> // std::sort
#include <iostream>

#include "state.hpp"

#include <ccd/collision_penalty_diff.hpp>
#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <io/opt_results.hpp>
#include <io/read_scene.hpp>
#include <io/write_scene.hpp>

#include <logger.hpp>

#include <autodiff/finitediff.hpp>
#include <opt/barrier_constraint.hpp>
#include <opt/displacement_opt.hpp>
#include <opt/ncp_solver.hpp>
#include <opt/volume_constraint.hpp>

namespace ccd {
State::State()
    : detection_method(DetectionMethod::BRUTE_FORCE)
    , volume_epsilon(1E-3)
    , output_dir(DATA_OUTPUT_DIR)
    , constraint_function(ConstraintType::VOLUME)
    , recompute_collision_set(false)
    , canvas_width(10)
    , canvas_height(10)
    , current_ev_impact(-1)
    , current_volume(-1)
    , grad_scaling(1.0f)
    , use_opt_gradient(false)
    , current_opt_time(0.0)
    , current_opt_iteration(-1)
{
}

// -----------------------------------------------------------------------------
// CRUD Scene
// -----------------------------------------------------------------------------

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);
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
    reset_impacts();

    current_volume = -1;
    current_ev_impact = -1;
    current_time = 0.0;
    current_opt_time = 0.0;
    selected_displacements.clear();
    selected_points.clear();

    opt_results.x.resizeLike(displacements);
    opt_results.x.setZero();
    opt_results.minf = 2e19;
    opt_results.success = false;
    opt_results.finished = false;

    opt_problem.fixed_dof.resize(displacements.size());
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

// -----------------------------------------------------------------------------
// CCD
// -----------------------------------------------------------------------------
void State::reset_impacts()
{
    volumes.resize(edges.rows());
    volumes.setZero();

    edge_impact_map.resize(edges.rows());
    edge_impact_map.setConstant(-1);

    ev_impacts.clear();
    ee_impacts.clear();

    volume_grad.resize(vertices.size(), edges.rows());
    volume_grad.setZero();
}

void State::run_ccd_pipeline()
{
    auto& collision_constraint = getCollisionConstraint();
    collision_constraint.initialize(vertices, edges, displacements);
    collision_constraint.compute_constraints(
        vertices, edges, displacements, volumes);
    collision_constraint.compute_constraints_jacobian(
        vertices, edges, displacements, volume_grad);

    current_volume = 0;
    current_ev_impact = ev_impacts.size() > 0 ? 0 : -1;

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

void State::reset_optimization_problem()
{

    volume_constraint.volume_epsilon = volume_epsilon;
    opt_problem.initialize(
        vertices, edges, displacements, getCollisionConstraint());
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

void State::optimize_displacements(const std::string filename)
{
    reset_optimization_problem();
    reset_results();

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
              u.resize(x.rows() / 2, 2);
              u_history.push_back(u);
              f_history.push_back(obj_value);

              Eigen::VectorXd gx = opt_problem.eval_g(x);
              g_history.push_back(gx);
              gsum_history.push_back(gx.sum());
              jac_g_history.push_back(opt_problem.eval_jac_g(x));

              it_x.push_back(x);
              it_lambda.push_back(dual);
              it_gamma.push_back(gamma);
          };

    // 5. run optimization
    if (solver_settings.method == ccd::opt::NCP) {
        opt_results = ncp_displ_solver.solve(opt_problem);
        opt_results.x.resize(U0.rows(), 2);
    } else if (solver_settings.method == ccd::opt::IPOPT) {
        ipopt_solver.max_iterations = solver_settings.max_iter;
        ipopt_solver.tolerance = solver_settings.relative_tolerance;
        ipopt_solver.print_level = solver_settings.verbosity;

        opt_results = ipopt_solver.solve(opt_problem);
        opt_results.x.resize(U0.rows(), 2);

    } else {
        opt_results = ccd::opt::displacement_optimization(
            opt_problem, U0, solver_settings);
    }

    if (solver_settings.verbosity > 0) {
        log_optimization_steps(filename, it_x, it_lambda, it_gamma);
    }

    opt_results.method = solver_settings.method;
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

void State::log_optimization_steps(const std::string filename,
    std::vector<Eigen::VectorXd>& it_x, std::vector<Eigen::VectorXd>& it_lambda,
    std::vector<double>& it_gamma)
{

    std::ofstream o;
    std::string sep = ",";
    if (filename.empty()) {
        o.copyfmt(std::cout);
        o.clear(std::cout.rdstate());
        o.basic_ios<char>::rdbuf(std::cout.rdbuf());
        sep = "\t";
    } else {
        o.open(filename);
    }

    // LOGGING for testing only --- -should move somewhere else
    const Eigen::IOFormat CommaFmt(
        Eigen::StreamPrecision, Eigen::DontAlignCols, sep, sep, "", "", "", "");

    // print table header
    o << ccd::opt::OptimizationMethodNames[solver_settings.method] << sep;
    o << ccd::opt::LCPSolverNames[solver_settings.lcp_solver] << std::endl;
    o << "it";
    for (int i = 0; i < opt_problem.num_vars; i++) {
        o << sep << "x_" << i;
    }
    for (int i = 0; i < opt_problem.num_constraints; i++) {
        o << sep << "lambda_" << i;
    }
    o << sep << "gamma" << sep << "f(x)";
    for (int i = 0; i < opt_problem.num_constraints; i++) {
        o << sep << "gx_" << i;
    }
    o << sep << "||jac_g(x)||" << std::endl;

    // print initial value
    Eigen::MatrixXd d = displacements;
    d.resize(d.size(), 1);
    o << 0 << sep;
    o << std::scientific << std::setprecision(2) << d.format(CommaFmt) << sep;
    o << Eigen::VectorXd::Zero(opt_problem.num_constraints).format(CommaFmt)
      << sep;
    o << 0.0 << sep;
    o << opt_problem.eval_f(d) << sep;
    o << opt_problem.eval_g(d).format(CommaFmt) << sep;

    Eigen::MatrixXd jac_gx = opt_problem.eval_jac_g(d);
    o << jac_gx.squaredNorm() << std::endl << std::endl;

    // print steps
    for (uint i = 0; i < it_x.size(); i++) {
        o << i + 1 << sep;
        o << it_x[i].format(CommaFmt) << sep;
        o << it_lambda[i].format(CommaFmt) << sep;
        o << it_gamma[i] << sep;
        o << opt_problem.eval_f(it_x[i]) << sep;
        o << opt_problem.eval_g(it_x[i]).format(CommaFmt) << sep;
        jac_gx = opt_problem.eval_jac_g(it_x[i]);
        o << jac_gx.squaredNorm() << std::endl;
    }
    spdlog::debug("Optimization Log saved to file://{}", filename);
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
    if (current_opt_iteration >= 0 && u_history.size() > 0) {
        current_opt_iteration %= u_history.size();
        displ = u_history[size_t(current_opt_iteration)];
    }
    return displ;
}

Eigen::MatrixX2d State::get_opt_volume_grad()
{
    if (current_volume < 0 || current_opt_iteration < 0
        || jac_g_history.size() == 0) {
        return Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    }
    auto& grad_volume = jac_g_history[size_t(current_opt_iteration)];
    if (grad_volume.size() == 0) {
        return Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    }
    current_volume %= grad_volume.rows();
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
    auto fun = opt_results.minf;
    if (current_opt_iteration >= 0 && f_history.size() > 0) {
        current_opt_iteration %= f_history.size();
        fun = f_history[size_t(current_opt_iteration)];
    }
    return fun;
}
} // namespace ccd
