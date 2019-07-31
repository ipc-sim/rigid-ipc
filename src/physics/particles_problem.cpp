#include <physics/particles_problem.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <io/serialize_json.hpp>
#include <opt/constraint_factory.hpp>
#include <physics/mass_matrix.hpp>
#include <utils/flatten.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace physics {

    ParticlesDisplProblem::ParticlesDisplProblem()
        : ParticlesDisplProblem("particle_problem")
    {
    }

    ParticlesDisplProblem::ParticlesDisplProblem(const std::string& name)
        : SimulationProblem(name)
        , constraint_ptr(nullptr)
        , intermediate_callback(nullptr)
        , use_mass_matrix(true)
        , collision_eps(2.0)
        , update_constraint_set(true)
        , is_linesearch_active(false)

    {
    }

    void ParticlesDisplProblem::init(const nlohmann::json& params)
    {

        constraint_ptr = opt::ConstraintFactory::factory().get_constraint(
            params["constraint"]);

        use_mass_matrix = params["use_mass_matrix"].get<bool>();
        update_constraint_set = params["update_constraint_set"].get<bool>();
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
        if (use_mass_matrix) {
            Eigen::VectorXd vertex_masses;
            physics::mass_vector(vertices_, edges_, vertex_masses);
            mass_matrix
                = Eigen::SparseDiagonal<double>(vertex_masses.replicate(2, 1));
        } else {
            mass_matrix = Eigen::SparseMatrix<double>(
                int(vertices_.size()), int(vertices_.size()));
            mass_matrix.setIdentity();
        }
        inv_mass_matrix = mass_matrix.cwiseInverse();

        // Collisions Groups
        // -------------------------------------------------
        // each particle is its own group
        group_ids_ = Eigen::VectorXi::LinSpaced(
            vertices_.rows(), 0, int(vertices_.rows()));

        // Simulation State
        // -------------------------------------------------
        Fcollision.resizeLike(vertices_);
        Fcollision.setZero();
        vertices_prev_ = vertices_;
        update_constraint();
    }
    nlohmann::json ParticlesDisplProblem::settings() const
    {
        nlohmann::json json;
        json["constraint"] = constraint_ptr->name();
        json["use_mass_matrix"] = use_mass_matrix;
        json["update_constraint_set"] = update_constraint_set;
        json["collision_eps"] = collision_eps;
        json["gravity"] = io::to_json(gravity_);
        return json;
    }

    bool ParticlesDisplProblem::simulation_step(const double time_step)
    {
        assert(constraint_ptr != nullptr);

        vertices_prev_ = vertices_;
        vertices_ = vertices_t1 = vertices_next(time_step);
        velocities_ = (vertices_ - vertices_prev_) / time_step;

        Fcollision.setZero();

        // check for collisions
        Eigen::MatrixXd& q0 = vertices_prev_;
        Eigen::MatrixXd& q1 = vertices_;

        return detect_collisions(q0, q1, CollisionCheck::CONSERVATIVE);
    }

    void ParticlesDisplProblem::update_constraint()
    {
        vertices_t0 = vec_vertices_t0 = vertices_prev_;
        vertices_t1 = vec_vertices_t1 = vertices_;

        flatten(vec_vertices_t0);
        flatten(vec_vertices_t1);

        constraint_ptr->initialize(
            vertices_t0, edges_, group_ids_, vertices_t1 - vertices_t0);

        // base problem initial solution
        x0 = vec_vertices_t0; // start from collision free state
        num_vars = int(x0.size());
    }

    bool ParticlesDisplProblem::take_step(
        const Eigen::VectorXd& x, const double time_step)
    {
        spdlog::trace(
            "particles_problem step=take_step xi={}", log::fmt_eigen(x));

        assert(constraint_ptr != nullptr);
        const Eigen::VectorXd& q1_collision = vec_vertices_t1;

        // update collision forces q* = q_collision + Fcollision
        //           Fc = (q* - q) M / (dt^2)
        Fcollision = mass_matrix * (x - q1_collision) / (time_step * time_step);
        ccd::unflatten(Fcollision, 2);

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

    /// ------------------------------------------------------------------------------
    /// GET state of system
    /// ------------------------------------------------------------------------------

    bool ParticlesDisplProblem::detect_collisions(const Eigen::MatrixXd& q0,
        const Eigen::MatrixXd& q1,
        const CollisionCheck check_type) const
    {
        assert(q0.cols() == 2);
        assert(q1.cols() == 2);

        EdgeVertexImpacts ev_impacts;
        double scale
            = check_type == CollisionCheck::EXACT ? 1.0 : (1.0 + collision_eps);

        ccd::detect_edge_vertex_collisions(q0, (q1 - q0) * scale, edges_,
            group_ids_, ev_impacts, constraint_ptr->detection_method,
            /*reset_impacts=*/true);
        return ev_impacts.size() > 0;
    }

    Eigen::MatrixXd ParticlesDisplProblem::vertices_next(
        const double time_step) const
    {
        Eigen::MatrixXd x = vertices_;
        x += time_step * velocities_; // momentum
        x.rowwise()
            += time_step * time_step * gravity_.transpose(); // body-forces
        x = (is_particle_dof_fixed).select(vertices_, x); // reset fixed nodes
        return x;
    }

    Eigen::MatrixXd ParticlesDisplProblem::velocities(
        const bool as_delta, const double time_step) const
    {
        if (as_delta) {
            return velocities_ * time_step;
        } else
            return velocities_;
    }
    Eigen::MatrixXd ParticlesDisplProblem::collision_force(
        const bool as_delta, const double time_step) const
    {
        if (as_delta) {
            Eigen::MatrixXd Fc_delta = Fcollision;
            ccd::flatten(Fc_delta);
            Fc_delta = inv_mass_matrix * Fc_delta * (time_step * time_step);
            ccd::unflatten(Fc_delta, 2);
            return Fc_delta;
        }
        return Fcollision;
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional
    ////////////////////////////////////////////////////////////////////////////

    double ParticlesDisplProblem::eval_f(const Eigen::VectorXd& q1)
    {
        Eigen::VectorXd diff = q1 - vec_vertices_t1;
        return 0.5 * diff.transpose() * mass_matrix * diff;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_grad_f(
        const Eigen::VectorXd& q1)
    {
        return mass_matrix * (q1 - vec_vertices_t1);
    }

    Eigen::SparseMatrix<double> ParticlesDisplProblem::eval_hessian_f(
        const Eigen::VectorXd& /*x*/)
    {
        return mass_matrix;
    }

    void ParticlesDisplProblem::eval_f_and_fdiff(const Eigen::VectorXd& x,
        double& f_uk,
        Eigen::VectorXd& f_uk_grad,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        f_uk = eval_f(x);
        f_uk_grad = eval_grad_f(x);
        f_uk_hessian = eval_hessian_f(x);
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Constraints
    ////////////////////////////////////////////////////////////////////////////
    Eigen::MatrixXd ParticlesDisplProblem::update_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = q1 - vec_vertices_t0;
        ccd::unflatten(Uk, 2);

        if (!is_linesearch_active && update_constraint_set) {
            constraint_ptr->detectCollisions(Uk);
        }
        return Uk;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);

        Eigen::VectorXd g_uk;
        constraint_ptr->compute_constraints(Uk, g_uk);
        return g_uk;
    };

    void ParticlesDisplProblem::eval_g(const Eigen::VectorXd& q1,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        constraint_ptr->compute_constraints(
            Uk, g_uk, g_uk_jacobian, g_uk_active);

#ifdef WITH_DERIVATIVE_CHECK
        double diff;
        if (!compare_jac_g_approx(q1, g_uk_jacobian.toDense(), diff)) {
            spdlog::error(
                "Derivative Check Failed in `eval_jac_g` diff_norm={} jac_norm={} g_sum={}",
                diff, g_uk_jacobian.norm(), g_uk.sum());
        }
#endif
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_jac_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        Eigen::MatrixXd jac_gx;
        constraint_ptr->compute_constraints_jacobian(Uk, jac_gx);

#ifdef WITH_DERIVATIVE_CHECK
        double diff;
        if (!compare_jac_g_approx(q1, jac_gx, diff)) {
            spdlog::error(
                "Derivative Check Failed in `eval_jac_g` diff_norm={}", diff);
        }
#endif
        return jac_gx;
    };

    void ParticlesDisplProblem::eval_jac_g(
        const Eigen::VectorXd& q1, Eigen::SparseMatrix<double>& jac_gx)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        constraint_ptr->compute_constraints_jacobian(Uk, jac_gx);

#ifdef WITH_DERIVATIVE_CHECK
        double diff;
        if (!compare_jac_g_approx(q1, jac_gx.toDense(), diff)) {
            spdlog::error(
                "Derivative Check Failed in `eval_jac_g` diff_norm={} jac_norm={}",
                diff, jac_gx.norm());
        }
#endif
    };

    std::vector<Eigen::SparseMatrix<double>>
    ParticlesDisplProblem::eval_hessian_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        constraint_ptr->compute_constraints_hessian(Uk, hess_gx);
        return hess_gx;
    }

    void ParticlesDisplProblem::eval_g_and_gdiff(const Eigen::VectorXd& q1,
        Eigen::VectorXd& g_uk,
        Eigen::MatrixXd& g_uk_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian)
    {
        Eigen::MatrixXd Uk = update_g(q1);

        constraint_ptr->compute_constraints_and_derivatives(
            Uk, g_uk, g_uk_jacobian, g_uk_hessian);
    }

    ////////////////////////////////////////////////////////////////////////////

    bool ParticlesDisplProblem::eval_intermediate_callback(
        const Eigen::VectorXd& x)
    {
        if (intermediate_callback != nullptr) {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);
            return intermediate_callback(x, Uk);
        }

        return true;
    }

    void ParticlesDisplProblem::enable_line_search_mode(
        const Eigen::VectorXd& max_x)
    {
        update_g(max_x);
        is_linesearch_active = true;
    }

    void ParticlesDisplProblem::disable_line_search_mode()
    {
        is_linesearch_active = false;
    }

} // namespace physics
} // namespace ccd
