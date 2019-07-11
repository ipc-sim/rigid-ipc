#include "rigid_body_problem.hpp"

#include <iostream>

#include <utils/flatten.hpp>
#include <utils/invalid_param_error.hpp>
#include <utils/tensor.hpp>

#include <io/read_rb_scene.hpp>
#include <opt/constraint_factory.hpp>

#include <logger.hpp>
namespace ccd {

namespace physics {

    RigidBodyProblem::RigidBodyProblem()
        : SimulationProblem("RigidBody")
        , use_chain_functional(false)
        , update_constraint_set(true)
    {
    }

    bool RigidBodyProblem::validate_params(const nlohmann::json& in_params)
    {
        std::vector<std::string> params = { { "rigid_bodies", "constraint",
            "use_chain_functional", "update_constraint_set" } };

        bool all_valid = true;
        for (auto& p : in_params.items()) {
            if (std::find(params.begin(), params.end(), p.key())
                == params.end()) {
                all_valid = false;
                spdlog::error("key {} not used", p.key());
                break;
            }
        }
        return all_valid;
    }

    void RigidBodyProblem::init(const nlohmann::json& params)
    {
        if (!validate_params(params)) {
            throw InvalidParameterError();
        }

        std::vector<physics::RigidBody> rbs;
        io::read_rb_scene(params, rbs);
        m_assembler.init(rbs);

        m_constraint_ptr = opt::ConstraintFactory::factory().get_constraint(
            params["constraint"]);

        // set parameters
        use_chain_functional = params["use_chain_functional"].get<bool>();
        update_constraint_set = params["update_constraint_set"].get<bool>();

        update_constraint();
    }

    void RigidBodyProblem::init(
        const std::vector<RigidBody> rbs, const std::string& constraint)
    {
        m_assembler.init(rbs);
        m_constraint_ptr
            = opt::ConstraintFactory::factory().get_constraint(constraint);

        update_constraint();
    }

    bool RigidBodyProblem::simulation_step(const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        for (auto& rb : m_assembler.m_rbs) {
            rb.position_prev = rb.position;
            rb.position += time_step * rb.velocity;
        }
        Eigen::MatrixXd q0 = m_assembler.world_vertices_t0();
        Eigen::MatrixXd q1 = m_assembler.world_vertices_t1();

        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(q0, q1 - q0, m_assembler.m_edges,
            ev_impacts, m_constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    bool RigidBodyProblem::take_step(
        const Eigen::VectorXd& rb_positions, const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        m_assembler.set_rb_positions(rb_positions);
        for (auto& rb : m_assembler.m_rbs) {
            rb.velocity = (rb.position - rb.position_prev) / time_step;
        }
        // check for collisions
        Eigen::MatrixXd q0 = m_assembler.world_vertices_t0();
        Eigen::MatrixXd q1 = m_assembler.world_vertices_t1();

        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(q0, q1 - q0, m_assembler.m_edges,
            ev_impacts, m_constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    void RigidBodyProblem::update_constraint()
    {
        m_q0 = m_assembler.world_vertices_t0();
        m_q1 = m_assembler.world_vertices_t1();

        m_sigma1 = m_assembler.rb_positions_t1();

        // base problem initial solution
        x0 = m_assembler.rb_positions_t0(); // start from collision free state
        num_vars = int(x0.size());

        // intended displacements
        Eigen::MatrixXd q1 = m_assembler.world_vertices_t1();
        m_constraint_ptr->initialize(m_q0, m_assembler.m_edges, q1 - m_q0);
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional
    ////////////////////////////////////////////////////////////////////////////
    double RigidBodyProblem::eval_f(const Eigen::VectorXd& sigma)
    {
        if (use_chain_functional) {
            return eval_f_chain(sigma);
        } else {
            Eigen::VectorXd diff = sigma - m_sigma1;
            return 0.5 * diff.transpose() * m_assembler.m_rb_mass_matrix * diff;
        }
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd grad_f;
        if (use_chain_functional) {
            grad_f = eval_grad_f_chain(sigma);
        } else {
            grad_f = m_assembler.m_rb_mass_matrix * (sigma - m_sigma1);
        }

#ifdef WITH_DERIVATIVE_CHECK
        assert(compare_grad_f_approx(sigma, grad_f));
#endif
        return grad_f;
    }

    Eigen::SparseMatrix<double> RigidBodyProblem::eval_hessian_f(
        const Eigen::VectorXd& sigma)
    {
        Eigen::SparseMatrix<double> hessian_f;
        if (use_chain_functional) {
            hessian_f = eval_hessian_f_chain(sigma);
        } else {
            hessian_f = m_assembler.m_rb_mass_matrix;
        }
#ifdef WITH_DERIVATIVE_CHECK
        assert(compare_hessian_f_approx(sigma, hessian_f));
#endif
        return hessian_f;
    }

    void RigidBodyProblem::eval_f_and_fdiff(const Eigen::VectorXd& sigma,
        double& f_uk,
        Eigen::VectorXd& f_uk_grad,
        Eigen::SparseMatrix<double>& f_uk_hessian)
    {
        // TODO: make this call efficient
        f_uk = eval_f(sigma);
        f_uk_grad = eval_grad_f(sigma);
        f_uk_hessian = eval_hessian_f(sigma);
#ifdef WITH_DERIVATIVE_CHECK
        assert(compare_grad_f_approx(sigma, f_uk_grad));
        assert(compare_hessian_f_approx(sigma, f_uk_hessian));
#endif
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional with CHAIN RULE
    ////////////////////////////////////////////////////////////////////////////
    double RigidBodyProblem::eval_f_chain(const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd qk = m_assembler.world_vertices(sigma);
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - m_q1));

        return 0.5
            * (diff.transpose() * m_assembler.m_mass_matrix * diff).sum();
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f_chain(
        const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd qk = m_assembler.world_vertices(sigma);
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - m_q1));

        Eigen::VectorXd grad_qk = m_assembler.m_mass_matrix * diff;

        Eigen::SparseMatrix<double> jac_qk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_qk_sigma);

        return grad_qk.transpose() * jac_qk_sigma;
    }

    Eigen::SparseMatrix<double> RigidBodyProblem::eval_hessian_f_chain(
        const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd qk = m_assembler.world_vertices(sigma);
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - m_q1));

        Eigen::SparseMatrix<double> jac_qk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_qk_sigma);
        std::vector<Eigen::SparseMatrix<double>> hess_qk_sigma;
        m_assembler.world_vertices_hessian(sigma, hess_qk_sigma);

        Eigen::VectorXd grad_qk = m_assembler.m_mass_matrix * diff;
        Eigen::SparseMatrix<double> hess_qk = m_assembler.m_mass_matrix;

        return jac_qk_sigma.transpose() * hess_qk * jac_qk_sigma
            + tensor::multiply(grad_qk, hess_qk_sigma);
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Constraints
    ////////////////////////////////////////////////////////////////////////////
    Eigen::MatrixXd RigidBodyProblem::update_g(const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd m_xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = m_xk - m_q0;

        if (update_constraint_set) {
            m_constraint_ptr->detectCollisions(uk);
        }
        return uk;
    }

    Eigen::VectorXd RigidBodyProblem::eval_g(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd g_uk;
        Eigen::MatrixXd uk = update_g(sigma);
        m_constraint_ptr->compute_constraints(uk, g_uk);
        return g_uk;
    }

    Eigen::MatrixXd RigidBodyProblem::eval_jac_g(const Eigen::VectorXd& sigma)
    {

        Eigen::MatrixXd uk = update_g(sigma);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);

        Eigen::MatrixXd jac_g_uk;
        m_constraint_ptr->compute_constraints_jacobian(uk, jac_g_uk);

        Eigen::MatrixXd jac = jac_g_uk * jac_xk_sigma;

#ifdef WITH_DERIVATIVE_CHECK
        assert(compare_jac_g_approx(sigma, jac));
#endif
        return jac;
    }

    /// @brief: util function to assemble hessian from partial derivatives
    void assemble_hessian(const Eigen::SparseMatrix<double>& jac_xk_sigma,
        const std::vector<Eigen::SparseMatrix<double>>& hess_xk_sigma,
        const Eigen::MatrixXd& jac_g_uk,
        const std::vector<Eigen::SparseMatrix<double>>& hessian_g_uk,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        size_t num_constraints = size_t(jac_g_uk.rows());
        assert(hessian_g_uk.size() == num_constraints);
        gx_hessian.clear();
        gx_hessian.reserve(num_constraints);

        for (unsigned long i = 0; i < num_constraints; i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(
                jac_xk_sigma.transpose() * hessian_g_uk[i] * jac_xk_sigma
                + tensor::multiply(jac_g_uk.row(long(i)), hess_xk_sigma));
        }
    }

    std::vector<Eigen::SparseMatrix<double>> RigidBodyProblem::eval_hessian_g(
        const Eigen::VectorXd& sigma)
    {

        Eigen::MatrixXd uk = update_g(sigma);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        std::vector<Eigen::SparseMatrix<double>> hess_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);
        m_assembler.world_vertices_hessian(sigma, hess_xk_sigma);

        Eigen::MatrixXd jac_g_uk;
        std::vector<Eigen::SparseMatrix<double>> hessian_g_uk;
        m_constraint_ptr->compute_constraints_jacobian(uk, jac_g_uk);
        m_constraint_ptr->compute_constraints_hessian(uk, hessian_g_uk);

        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        assemble_hessian(
            jac_xk_sigma, hess_xk_sigma, jac_g_uk, hessian_g_uk, gx_hessian);
        return gx_hessian;
    }

    void RigidBodyProblem::eval_g_and_gdiff(const Eigen::VectorXd& sigma,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& gx_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        Eigen::MatrixXd uk = update_g(sigma);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        std::vector<Eigen::SparseMatrix<double>> hess_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);
        m_assembler.world_vertices_hessian(sigma, hess_xk_sigma);

        Eigen::MatrixXd jac_g_uk;
        std::vector<Eigen::SparseMatrix<double>> hessian_g_uk;
        m_constraint_ptr->compute_constraints_and_derivatives(
            uk, gx, jac_g_uk, hessian_g_uk);

        gx_jacobian = jac_g_uk * jac_xk_sigma;

        assemble_hessian(
            jac_xk_sigma, hess_xk_sigma, jac_g_uk, hessian_g_uk, gx_hessian);

#ifdef WITH_DERIVATIVE_CHECK
        assert(compare_jac_g_approx(sigma, gx_jacobian));
#endif
    }

} // namespace physics
} // namespace ccd
