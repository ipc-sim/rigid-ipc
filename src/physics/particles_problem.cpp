#include <physics/particles_problem.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/serialize_json.hpp>
#include <physics/mass_matrix.hpp>
#include <utils/flatten.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace physics {

    ParticlesProblem::ParticlesProblem()
        : ParticlesProblem("particle_problem")
    {
    }

    ParticlesProblem::ParticlesProblem(const std::string& name)
        : collision_eps(2.0)
        , name_(name)

    {
    }

    void ParticlesProblem::settings(const nlohmann::json& params)
    {

        collision_eps = params["collision_eps"].get<double>();

        io::from_json(params["vertices"], vertices_);
        io::from_json(params["edges"], edges_);
        io::from_json(params["velocities"], velocities_);
        assert(vertices_.rows() == velocities_.rows());

        io::from_json(params["gravity"], gravity_);
        assert(gravity_.rows() == 2);

        // some degrees of freedom are actually fixed
        // -------------------------------------------------
        is_particle_dof_fixed.resize(vertices_.rows(), 2);
        is_particle_dof_fixed.setConstant(false);
        for (auto i : params["x_fixed"].get<std::vector<int>>()) {
            is_particle_dof_fixed(i, 0) = true;
        }
        for (auto i : params["y_fixed"].get<std::vector<int>>()) {
            is_particle_dof_fixed(i, 1) = true;
        }
        is_dof_fixed_ = flat<bool>(is_particle_dof_fixed);

        // mass matrix
        // -------------------------------------------------

        Eigen::VectorXd vertex_masses;
        physics::mass_vector(vertices_, edges_, vertex_masses);
        mass_matrix
            = Eigen::SparseDiagonal<double>(vertex_masses.replicate(2, 1));


        // Collisions Groups
        // -------------------------------------------------
        // each particle is its own group
        group_ids_ = Eigen::VectorXi::LinSpaced(
            vertices_.rows(), 0, int(vertices_.rows()));

        // Simulation State
        // -------------------------------------------------
        vertices_prev_ = vertices_;
        update_constraint();
    }
    nlohmann::json ParticlesProblem::settings() const
    {
        nlohmann::json json;
        json["collision_eps"] = collision_eps;
        json["gravity"] = io::to_json(gravity_);
        return json;
    }

    nlohmann::json ParticlesProblem::state() const
    {
        nlohmann::json json;
        json["vertices"] = io::to_json(vertices_);
        json["velocities"] = io::to_json(velocities_);
        return json;
    }
    void ParticlesProblem::state(const nlohmann::json& params)
    {
        nlohmann::json json;

        io::from_json(params["vertices"], vertices_);
        io::from_json(params["velocities"], velocities_);
    }

    bool ParticlesProblem::simulation_step(const double time_step)
    {
        vertices_prev_ = vertices_;
        vertices_ = vertices_t1 = vertices_next(time_step);
        velocities_ = (vertices_ - vertices_prev_) / time_step;

        // check for collisions
        Eigen::MatrixXd& q0 = vertices_prev_;
        Eigen::MatrixXd& q1 = vertices_;

        return detect_collisions(q0, q1, CollisionCheck::CONSERVATIVE);
    }

    void ParticlesProblem::update_constraint()
    {
        vertices_t0 = vec_vertices_t0 = vertices_prev_;
        vertices_t1 = vec_vertices_t1 = vertices_;

        flatten(vec_vertices_t0);
        flatten(vec_vertices_t1);

        constraint().initialize(
            vertices_t0, edges_, group_ids_, vertices_t1 - vertices_t0);

        // base problem initial solution
        x0 = vec_vertices_t0; // start from collision free state
        num_vars_ = int(x0.size());
    }

    opt::OptimizationResults ParticlesProblem::solve_constraints()
    {
        return solver().solve();
    }

    void ParticlesProblem::init_solve() { return solver().init_solve(); }
    opt::OptimizationResults ParticlesProblem::step_solve()
    {
        return solver().step_solve();
    }

    bool ParticlesProblem::take_step(
        const Eigen::VectorXd& x, const double time_step)
    {
        spdlog::trace(
            "particles_problem step=take_step xi={}", log::fmt_eigen(x));

        // update final position
        vertices_ = x;

        // update velocities
        ccd::unflatten(vertices_, 2);
        velocities_ = (vertices_ - vertices_prev_) / time_step;

        // check for collisions
        Eigen::MatrixXd& q0 = vertices_prev_;
        Eigen::MatrixXd& q1 = vertices_;
        return detect_collisions(q0, q1, CollisionCheck::EXACT);
    }

    bool ParticlesProblem::detect_collisions(const Eigen::MatrixXd& q0,
        const Eigen::MatrixXd& q1,
        const CollisionCheck check_type)
    {
        assert(q0.cols() == 2);
        assert(q1.cols() == 2);

        EdgeVertexImpacts ev_impacts;
        double scale
            = check_type == CollisionCheck::EXACT ? 1.0 : (1.0 + collision_eps);

        ccd::detect_edge_vertex_collisions(q0, (q1 - q0) * scale, edges_,
            group_ids_, ev_impacts, constraint().detection_method);
        return ev_impacts.size() > 0;
    }

    Eigen::MatrixXd ParticlesProblem::vertices_next(
        const double time_step) const
    {
        Eigen::MatrixXd x = vertices_;
        x += time_step * velocities_; // momentum
        x.rowwise()
            += time_step * time_step * gravity_.transpose(); // body-forces
        x = (is_particle_dof_fixed).select(vertices_, x); // reset fixed nodes
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional
    ////////////////////////////////////////////////////////////////////////////

    double ParticlesProblem::eval_f(const Eigen::VectorXd& q1)
    {
        Eigen::VectorXd diff = q1 - vec_vertices_t1;
        return 0.5 * diff.transpose() * mass_matrix * diff;
    }

    Eigen::VectorXd ParticlesProblem::eval_grad_f(const Eigen::VectorXd& q1)
    {
        return mass_matrix * (q1 - vec_vertices_t1);
    }

    Eigen::SparseMatrix<double> ParticlesProblem::eval_hessian_f(
        const Eigen::VectorXd& /*x*/)
    {
        return mass_matrix;
    }

} // namespace physics
} // namespace ccd
