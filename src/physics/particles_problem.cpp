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
        : SimulationProblem("ParticlesProblem")
        , constraint_ptr(nullptr)
        , intermediate_callback(nullptr)
        , use_mass_matrix(true)
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

        io::from_json(params["vertices"], vertices_);
        io::from_json(params["edges"], edges_);
        io::from_json(params["velocities"], velocity);
        assert(vertices_.rows() == velocity.rows());

        // some degrees of freedom are actually fixed
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

        // Simulation State
        // -------------------------------------------------
        vertices_prev.resize(vertices_.rows(), vertices_.cols());
        vertices_prev.setZero();

        update_constraint();
    }

    bool ParticlesDisplProblem::simulation_step(const double time_step)
    {
        assert(constraint_ptr != nullptr);

        vertices_prev = vertices_;
        vertices_ += time_step * velocity;

        // check for collisions
        Eigen::MatrixXd& q0 = vertices_prev;
        Eigen::MatrixXd& q1 = vertices_;
        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(q0, q1 - q0, edges_, ev_impacts,
            constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    bool ParticlesDisplProblem::take_step(
        const Eigen::VectorXd& x, const double time_step)
    {
        assert(constraint_ptr != nullptr);
        vertices_ = x;
        unflatten(vertices_, 2);

        velocity = (vertices_ - vertices_prev) / time_step;

        // check for collisions
        Eigen::MatrixXd& q0 = vertices_prev;
        Eigen::MatrixXd& q1 = vertices_;
        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(q0, q1 - q0, edges_, ev_impacts,
            constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    void ParticlesDisplProblem::update_constraint()
    {
        q0_ = vertices_prev;
        q1_ = vertices_;

        constraint_ptr->initialize(q0_, edges_, q1_ - q0_);

        flatten(q0_);
        flatten(q1_);

        // base problem initial solution
        x0 = q0_; // start from collision free state
        num_vars = int(x0.size());
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional
    ////////////////////////////////////////////////////////////////////////////

    double ParticlesDisplProblem::eval_f(const Eigen::VectorXd& q1)
    {
        Eigen::VectorXd diff = q1 - q1_;
        return 0.5 * diff.transpose() * mass_matrix * diff;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_grad_f(
        const Eigen::VectorXd& q1)
    {
        return mass_matrix * (q1 - q1_);
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
        Eigen::MatrixXd Uk = q1 - q0_;
        unflatten(Uk, 2);

        if (!is_linesearch_active && update_constraint_set) {
            constraint_ptr->detectCollisions(Uk);
        }
        return Uk;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);

        Eigen::VectorXd g_uk;
        PROFILE(constraint_ptr->compute_constraints(Uk, g_uk),
            ProfiledPoint::COMPUTING_CONSTRAINTS);
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
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_jac_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        Eigen::MatrixXd jac_gx;
        PROFILE(constraint_ptr->compute_constraints_jacobian(Uk, jac_gx),
            ProfiledPoint::COMPUTING_GRADIENT);

        return jac_gx;
    };

    void ParticlesDisplProblem::eval_jac_g(
        const Eigen::VectorXd& q1, Eigen::SparseMatrix<double>& jac_gx)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        constraint_ptr->compute_constraints_jacobian(Uk, jac_gx);
    };

    std::vector<Eigen::SparseMatrix<double>>
    ParticlesDisplProblem::eval_hessian_g(const Eigen::VectorXd& q1)
    {
        Eigen::MatrixXd Uk = update_g(q1);
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        PROFILE(constraint_ptr->compute_constraints_hessian(Uk, hess_gx);
                , ProfiledPoint::COMPUTING_HESSIAN);
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
