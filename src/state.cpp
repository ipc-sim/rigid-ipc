#include <algorithm> // std::sort
#include <iostream>

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
    : output_dir(DATA_OUTPUT_DIR)
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

    // Resize the fixed_dof
    Eigen::Array<bool, Eigen::Dynamic, 1> new_fixed_dof(displacements.size());
    new_fixed_dof.block(0, 0, lastid, 1)
        = opt_problem.fixed_dof.block(0, 0, lastid, 1);
    new_fixed_dof(lastid) = false;
    new_fixed_dof.block(lastid + 1, 0, lastid, 1)
        = opt_problem.fixed_dof.block(lastid, 0, lastid, 1);
    new_fixed_dof(new_fixed_dof.size() - 1) = false;
    opt_problem.fixed_dof = new_fixed_dof;
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);

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
    volumes.resize(0);
    volumes.setZero();

    volume_grad.resize(0, 0);
    volume_grad.setZero();
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
    f_history.push_back(opt_problem.eval_f(x));

    // Allways use Volume constraint
    volume_constraint.detecteCollisions(Uk);

    Eigen::VectorXd gx;
    Eigen::MatrixXd jac_gx;
    volume_constraint.compute_constraints(Uk, gx, jac_gx);
    g_history.push_back(gx);
    gsum_history.push_back(gx.sum());
    jac_g_history.push_back(jac_gx);

    // check if there is any reason to stop optimization
    return true;
}

void State::optimize_displacements(const std::string filename)
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
    Eigen::MatrixX2d empty_grad = Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    if (current_opt_iteration < 0 || jac_g_history.size() == 0
        || current_volume < 0) {
        return empty_grad;
    }

    auto& grad_volume = jac_g_history[size_t(current_opt_iteration)];
    if (grad_volume.size() == 0) {
        return empty_grad;
    }

    current_volume = std::max(-1, current_volume);
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
    auto fun = opt_results.minf;
    if (current_opt_iteration >= 0 && f_history.size() > 0) {
        current_opt_iteration %= f_history.size();
        fun = f_history[size_t(current_opt_iteration)];
    }
    return fun;
}
} // namespace ccd
