#include "rigid_body_problem.hpp"

#include <iostream>

#include <utils/flatten.hpp>
#include <utils/invalid_param_error.hpp>
#include <utils/tensor.hpp>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>

#include <opt/constraint_factory.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {

namespace physics {

    RigidBodyProblem::RigidBodyProblem()
        : RigidBodyProblem("rigid_body_problem")
    {
    }

    RigidBodyProblem::RigidBodyProblem(const std::string& name)
        : SimulationProblem(name)
        , use_chain_functional(false)
        , update_constraint_set(true)
        , gravity_(Eigen::Vector3d::Zero())
        , collision_eps(2)
    {
    }

    void RigidBodyProblem::init(const nlohmann::json& params)
    {
        std::vector<physics::RigidBody> rbs;
        io::read_rb_scene(params, rbs);
        m_assembler.init(rbs);

        m_constraint_ptr = opt::ConstraintFactory::factory().get_constraint(
            params["constraint"]);

        // set parameters
        use_chain_functional = params["use_chain_functional"].get<bool>();
        update_constraint_set = params["update_constraint_set"].get<bool>();
        collision_eps = params["collision_eps"].get<double>();
        io::from_json(params["gravity"], gravity_);
        assert(gravity_.rows() == 3);

        vertices_t0.resize(m_assembler.num_vertices(), 2);
        vertices_t0.setZero();

        vertices_q1.resize(m_assembler.num_vertices(), 2);
        vertices_q1.setZero();

        Fcollision.resize(m_assembler.num_vertices(), 2);
        Fcollision.setZero();
        update_constraint();
    }

    nlohmann::json RigidBodyProblem::settings() const
    {
        nlohmann::json json;
        json["constraint"] = m_constraint_ptr->name();
        json["use_chain_functional"] = use_chain_functional;
        json["update_constraint_set"] = update_constraint_set;
        json["collision_eps"] = collision_eps;
        json["gravity"] = io::to_json(gravity_);
        return json;
    }
    void RigidBodyProblem::init(
        const std::vector<RigidBody> rbs, const std::string& constraint)
    {
        m_assembler.init(rbs);
        m_constraint_ptr
            = opt::ConstraintFactory::factory().get_constraint(constraint);

        vertices_t0.resize(m_assembler.num_vertices(), 2);
        vertices_t0.setZero();

        vertices_q1.resize(m_assembler.num_vertices(), 2);
        vertices_q1.setZero();

        Fcollision.resize(m_assembler.num_vertices(), 2);
        Fcollision.setZero();
        update_constraint();
    }

    nlohmann::json RigidBodyProblem::state() const
    {
        nlohmann::json json;
        std::vector<nlohmann::json> rbs;
        for (auto& rb : m_assembler.m_rbs) {
            nlohmann::json jrb;
            jrb["position"] = io::to_json(Eigen::VectorXd(rb.position));
            jrb["velocity"] = io::to_json(Eigen::VectorXd(rb.velocity));
            rbs.push_back(jrb);
        }
        json["rigid_bodies"] = rbs;
        return json;
    }

    void RigidBodyProblem::state(const nlohmann::json& args)
    {
        nlohmann::json json;
        auto& rbs = args["rigid_bodies"];
        assert(rbs.size() == m_assembler.m_rbs.size());
        size_t i = 0;
        for (auto& jrb: args["rigid_bodies"] ) {
            Eigen::VectorXd velocity;
            io::from_json(jrb["velocity"], velocity);
            Eigen::VectorXd position;
            io::from_json(jrb["position"], position);

            m_assembler.m_rbs[i].position = position;
            m_assembler.m_rbs[i].velocity = velocity;
            i+=1;
        }

    }

    bool RigidBodyProblem::simulation_step(const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        for (auto& rb : m_assembler.m_rbs) {
            rb.position_prev = rb.position;
            rb.position = rb_position_next(rb, time_step);
            rb.velocity = (rb.position - rb.position_prev) / time_step;
        }

        Fcollision.setZero();
        vertices_q1.setZero();

        vertices_t0 = m_assembler.world_vertices_t0();
        vertices_q1 = m_assembler.world_vertices_t1();

        return detect_collisions(
            vertices_t0, vertices_q1, CollisionCheck::CONSERVATIVE);
    }

    bool RigidBodyProblem::take_step(
        const Eigen::VectorXd& rb_positions, const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        const Eigen::MatrixXd& q1_collisions = vertices_q1;

        // update final position
        m_assembler.set_rb_positions(rb_positions);
        Eigen::MatrixXd q1_new = m_assembler.world_vertices_t1();

        // update collision forces // q_new = q_collision + Fcollision
        Eigen::MatrixXd delta_c = (q1_new - q1_collisions);
        flatten(delta_c);
        Fcollision
            = m_assembler.m_mass_matrix * delta_c / (time_step * time_step);
        ccd::unflatten(Fcollision, 2);

        // update velocities
        for (auto& rb : m_assembler.m_rbs) {
            rb.velocity = (rb.position - rb.position_prev) / time_step;
        }

        return detect_collisions(vertices_t0, q1_new, CollisionCheck::EXACT);
    }
    void RigidBodyProblem::update_constraint()
    {
        vertices_t0 = m_assembler.world_vertices_t0();
        vertices_q1 = m_assembler.world_vertices_t1();

        rb_positions_t1 = m_assembler.rb_positions_t1();

        // base problem initial solution
        x0 = m_assembler.rb_positions_t0(); // start from collision free state
        num_vars = int(x0.size());

        m_constraint_ptr->initialize(vertices_t0, m_assembler.m_edges,
            m_assembler.m_vertex_to_body_map, vertices_q1 - vertices_t0);
    }

    /// -----------------------------------------------------------------

    Eigen::MatrixXd RigidBodyProblem::vertices_next(
        const double time_step) const
    {
        std::vector<Eigen::Vector3d> positions;
        positions.reserve(m_assembler.m_rbs.size());

        for (auto& rb : m_assembler.m_rbs) {
            Eigen::Vector3d pos_new;
            pos_new = rb_position_next(rb, time_step);
            positions.push_back(pos_new);
        }

        Eigen::MatrixXd q1 = m_assembler.world_vertices(positions);
        return q1;
    }

    Eigen::MatrixXd RigidBodyProblem::velocities(
        const bool as_delta, const double time_step) const
    {
        Eigen::MatrixXd vel = m_assembler.world_velocities();
        if (as_delta) {
            vel *= time_step;
        }
        return vel;
    }
    Eigen::MatrixXd RigidBodyProblem::collision_force(
        const bool as_delta, const double time_step) const
    {
        if (as_delta) {
            Eigen::MatrixXd Fc_delta = Fcollision;
            flatten(Fc_delta);
            Fc_delta = m_assembler.m_inv_mass_matrix * Fc_delta
                * (time_step * time_step);
            ccd::unflatten(Fc_delta, 2);
            return Fc_delta;
        } else {
            return Fcollision;
        }
    }

    Eigen::Vector3d RigidBodyProblem::rb_position_next(
        const RigidBody& rb, const double time_step) const
    {
        Eigen::Vector3d x = rb.position;
        x += time_step * rb.velocity;                 // momentum
        x += time_step * time_step * gravity_;        // body-forces
        x = (rb.is_dof_fixed).select(rb.position, x); // reset fixed nodes
        return x;
    }

    bool RigidBodyProblem::detect_collisions(const Eigen::MatrixXd& q0,
        const Eigen::MatrixXd& q1,
        const CollisionCheck check_type)
    {
        assert(q0.cols() == 2);
        assert(q1.cols() == 2);

        EdgeVertexImpacts ev_impacts;
        double scale
            = check_type == CollisionCheck::EXACT ? 1.0 : (1.0 + collision_eps);
        ccd::detect_edge_vertex_collisions(q0, (q1 - q0) * scale,
            m_assembler.m_edges, m_assembler.m_vertex_to_body_map, ev_impacts,
            m_constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    ////////////////////////////////////////////////////////////////////////////
    /// Functional
    ////////////////////////////////////////////////////////////////////////////
    double RigidBodyProblem::eval_f(const Eigen::VectorXd& sigma)
    {
        if (use_chain_functional) {
            return eval_f_chain(sigma);
        } else {
            Eigen::VectorXd diff = sigma - rb_positions_t1;
            return 0.5 * diff.transpose() * m_assembler.m_rb_mass_matrix * diff;
        }
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd grad_f;
        if (use_chain_functional) {
            grad_f = eval_grad_f_chain(sigma);
        } else {
            Eigen::VectorXd diff = sigma - rb_positions_t1;
            grad_f = m_assembler.m_rb_mass_matrix * diff;
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
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - vertices_q1));

        return 0.5
            * (diff.transpose() * m_assembler.m_mass_matrix * diff).sum();
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f_chain(
        const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd qk = m_assembler.world_vertices(sigma);
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - vertices_q1));

        Eigen::VectorXd grad_qk = m_assembler.m_mass_matrix * diff;

        Eigen::SparseMatrix<double> jac_qk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_qk_sigma);

        return grad_qk.transpose() * jac_qk_sigma;
    }

    Eigen::SparseMatrix<double> RigidBodyProblem::eval_hessian_f_chain(
        const Eigen::VectorXd& sigma)
    {
        Eigen::MatrixXd qk = m_assembler.world_vertices(sigma);
        Eigen::VectorXd diff = flat(Eigen::MatrixXd(qk - vertices_q1));

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
        Eigen::MatrixXd uk = m_xk - vertices_t0;

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
        bool aux = update_constraint_set;
        update_constraint_set = false;
        if (!compare_jac_g_approx(sigma, jac)) {
            spdlog::error(
                "rigid_body status=fail message='constraint finite-differences failed'");
        }
        update_constraint_set = aux;
#endif
        return jac;
    }

    void RigidBodyProblem::eval_jac_g(
        const Eigen::VectorXd& sigma, Eigen::SparseMatrix<double>& jac_gx)
    {
        Eigen::MatrixXd uk = update_g(sigma);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);

        m_constraint_ptr->compute_constraints_jacobian(uk, jac_gx);
        jac_gx = jac_gx * jac_xk_sigma;
    };

    void RigidBodyProblem::eval_g(const Eigen::VectorXd& sigma,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        Eigen::MatrixXd uk = update_g(sigma);

        Eigen::SparseMatrix<double> jac_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);

        m_constraint_ptr->compute_constraints(
            uk, g_uk, g_uk_jacobian, g_uk_active);

        g_uk_jacobian = g_uk_jacobian * jac_xk_sigma;
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
        NAMED_PROFILE_POINT("rigid_body_problem__update_g", UPDATE_G)
        NAMED_PROFILE_POINT(
            "rigid_body_problem__rigid_body_gradient", RIGID_BODY_GRADS)
        NAMED_PROFILE_POINT(
            "rigid_body_problem__rigid_body_hessian", RIGID_BODY_HESSIAN)
        NAMED_PROFILE_POINT(
            "rigid_body_problem__particles_gradients", PARTICLES_GRADS)
        NAMED_PROFILE_POINT(
            "rigid_body_problem__assemble_hessian", ASSEMBLE_HESSIAN)

        PROFILE_START(UPDATE_G)
        Eigen::MatrixXd uk = update_g(sigma);
        PROFILE_END(UPDATE_G)

        Eigen::SparseMatrix<double> jac_xk_sigma;
        std::vector<Eigen::SparseMatrix<double>> hess_xk_sigma;
        PROFILE_START(RIGID_BODY_GRADS)
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);
        PROFILE_END(RIGID_BODY_GRADS)
        PROFILE_START(RIGID_BODY_HESSIAN)
        m_assembler.world_vertices_hessian(sigma, hess_xk_sigma);
        PROFILE_END(RIGID_BODY_HESSIAN)

        Eigen::MatrixXd jac_g_uk;
        std::vector<Eigen::SparseMatrix<double>> hessian_g_uk;
        PROFILE_START(PARTICLES_GRADS)
        m_constraint_ptr->compute_constraints_and_derivatives(
            uk, gx, jac_g_uk, hessian_g_uk);
        PROFILE_END(PARTICLES_GRADS)

        gx_jacobian = jac_g_uk * jac_xk_sigma;
        PROFILE_START(ASSEMBLE_HESSIAN)
        assemble_hessian(
            jac_xk_sigma, hess_xk_sigma, jac_g_uk, hessian_g_uk, gx_hessian);
        PROFILE_END(ASSEMBLE_HESSIAN)

#ifdef WITH_DERIVATIVE_CHECK
        bool aux = update_constraint_set;
        update_constraint_set = false;
        if (!compare_jac_g_approx(sigma, gx_jacobian)) {
            spdlog::error(
                "rigid_body status=fail message='constraint finite-differences failed'");
        }
        update_constraint_set = aux;
#endif
    }

} // namespace physics
} // namespace ccd
