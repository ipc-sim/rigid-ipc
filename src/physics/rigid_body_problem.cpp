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
        , coefficient_restitution(0)
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
        coefficient_restitution
            = params["coefficient_restitution"].get<double>();
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
        json["coefficient_restitution"] = coefficient_restitution;
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
        Eigen::Vector2d p = Eigen::Vector2d::Zero();
        double L = 0;
        double T = 0;
        double G = 0;
        for (auto& rb : m_assembler.m_rbs) {
            nlohmann::json jrb;
            jrb["position"] = io::to_json(Eigen::VectorXd(rb.position));
            jrb["velocity"] = io::to_json(Eigen::VectorXd(rb.velocity));
            rbs.push_back(jrb);
            p += rb.mass * rb.velocity.head(2);
            L += rb.moment_of_inertia * rb.velocity(2);
            T += 1.0/2.0 * rb.mass * rb.velocity.head(2).transpose() * rb.velocity.head(2);
            T += 1.0/2.0 * rb.moment_of_inertia * rb.velocity(2) * rb.velocity(2);
            G += - rb.mass * gravity().transpose() * rb.position;

        }
        json["rigid_bodies"] = rbs;
        json["linear_momentum"] = io::to_json(Eigen::VectorXd(p));
        json["angular_momentum"] = L;
        json["kinetic_energy"] = T;
        json["potential_energy"] = G;
        return json;
    }

    void RigidBodyProblem::state(const nlohmann::json& args)
    {
        nlohmann::json json;
        auto& rbs = args["rigid_bodies"];
        assert(rbs.size() == m_assembler.m_rbs.size());
        size_t i = 0;
        for (auto& jrb : args["rigid_bodies"]) {
            Eigen::VectorXd velocity;
            io::from_json(jrb["velocity"], velocity);
            Eigen::VectorXd position;
            io::from_json(jrb["position"], position);

            m_assembler.m_rbs[i].position = position;
            m_assembler.m_rbs[i].velocity = velocity;
            i += 1;
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

        // update final position
        // -------------------------------------
        m_assembler.set_rb_positions(rb_positions);
        Eigen::MatrixXd q1 = m_assembler.world_vertices_t1();

        // [DEBUG!] update collision forces
        // -------------------------------------
        const Eigen::MatrixXd& q1_0 = vertices_q1;
        Eigen::MatrixXd delta_c = (q1 - q1_0);
        flatten(delta_c);
        Fcollision
            = m_assembler.m_mass_matrix * delta_c / (time_step * time_step);
        ccd::unflatten(Fcollision, 2);

        // update velocities
        // -------------------------------------
        solve_velocities();
        // TODO: check this two are the same !!!
//        if (coefficient_restitution > 0) {
//            solve_velocities();
//        } else {
//            for (auto& rb : m_assembler.m_rbs) {
//                rb.velocity = (rb.position - rb.position_prev) / time_step;
//            }
//        }

        return detect_collisions(vertices_t0, q1, CollisionCheck::EXACT);
    }

    void RigidBodyProblem::solve_velocities()
    {

        // precompute normal directions (since velocities will change i can't do
        // it after
        Eigen::MatrixXd normals(original_ev_impacts.size(), 2);

        for (int i = 0; i < normals.rows(); ++i) {
            const auto& ev_candidate = original_ev_impacts[size_t(i)];

            double toi = ev_candidate.time;

            const long edge_id = ev_candidate.edge_index;
            const long a_id = ev_candidate.vertex_index;
            const int b0_id = m_assembler.m_edges.coeff(edge_id, 0);
            const int b1_id = m_assembler.m_edges.coeff(edge_id, 1);

            Eigen::Vector2d n_toi;
            auto v_toi = vertices_t0 + (vertices_q1 - vertices_t0) * toi;
            Eigen::VectorXd e_toi
                = v_toi.row(b1_id) - v_toi.row(b0_id); // edge at toi
            n_toi << -e_toi(1), e_toi(0);              // 90deg rotation
            n_toi.normalize();

            // check correct diraction
            Eigen::Vector2d va = vertices_t0.row(a_id);
            Eigen::Vector2d vb = vertices_t0.row(b0_id);
            if ((va - vb).transpose() * n_toi <= 0.0) {
                n_toi *= -1;
            }
            normals.row(i) = n_toi.transpose();
        }

#ifndef NDEBUG
        double prev_toi = -1;
#endif
        for (int i = 0; i < normals.rows(); ++i) {
            const auto& ev_candidate = original_ev_impacts[size_t(i)];

            double toi = ev_candidate.time;
#ifndef NDEBUG
            assert(prev_toi < toi);
            prev_toi = toi;
#endif
            double alpha = ev_candidate.alpha;

            // global ids of the vertices
            const long edge_id = ev_candidate.edge_index;
            const long a_id = ev_candidate.vertex_index;
            const int b0_id = m_assembler.m_edges.coeff(edge_id, 0);
            const int b1_id = m_assembler.m_edges.coeff(edge_id, 1);

            const size_t body_A_id
                = size_t(m_assembler.m_vertex_to_body_map(a_id));
            const size_t body_B_id
                = size_t(m_assembler.m_vertex_to_body_map(b0_id));

            // local (rigid body) ids of the vertices
            const long r_A_id = a_id - m_assembler.m_body_vertex_id[body_A_id];
            const long r_B0_id
                = b0_id - m_assembler.m_body_vertex_id[body_B_id];
            const long r_B1_id
                = b1_id - m_assembler.m_body_vertex_id[body_B_id];

            auto& body_A = m_assembler.m_rbs[body_A_id];
            auto& body_B = m_assembler.m_rbs[body_B_id];

            // The velocities of the center of mass
            const Eigen::Vector2d& V_Aprev = body_A.velocity.head(2);
            const Eigen::Vector2d& V_Bprev = body_B.velocity.head(2);
            // The angular velocities
            const double& w_Aprev = body_A.velocity(2);
            const double& w_Bprev = body_B.velocity(2);
            // The masss
            const double inv_m_A
                = body_A.is_dof_fixed[0] || body_A.is_dof_fixed[1]
                ? 0.0
                : 1.0 / body_A.mass;
            const double inv_m_B
                = body_B.is_dof_fixed[0] || body_B.is_dof_fixed[1]
                ? 0.0
                : 1.0 / body_B.mass;
            // The moment of inertia
            const double inv_I_A
                = body_A.is_dof_fixed[2] ? 0 : 1.0 / body_A.moment_of_inertia;
            const double inv_I_B
                = body_B.is_dof_fixed[2] ? 0 : 1.0 / body_B.moment_of_inertia;

            // Vectors from the center of mass to the collision point
            // (90deg rotation counter clockwise)
            //
            //      (1) first get vertices position wrt rigid bodies
            const Eigen::Vector2d r0_A = body_A.vertices.row(r_A_id);
            const Eigen::Vector2d r0_B0
                = body_B.vertices.row(r_B0_id); // edge vertex 0
            const Eigen::Vector2d r0_B1
                = body_B.vertices.row(r_B1_id); // edge vertex 1
            const Eigen::Vector2d r0_B = r0_B0 + alpha * (r0_B1 - r0_B0);

            //      (2) and the angular displacement at time of collision
            const double theta_Atoi = body_A.position_prev(2)
                + toi * (body_A.position - body_A.position_prev)(2);
            const double theta_Btoi = body_B.position_prev(2)
                + toi * (body_B.position - body_B.position_prev)(2);
            //      (3) then the vectors are given by r = R(\theta_{t})*r_0
            const Eigen::Vector2d r_Aperp_toi
                = body_A.grad_theta(theta_Atoi) * r0_A;
            const Eigen::Vector2d r_Bperp_toi
                = body_A.grad_theta(theta_Btoi) * r0_B;

            // The collision point velocities BEFORE collision
            const Eigen::Vector2d v_Aprev = V_Aprev + w_Aprev * r_Aperp_toi;
            const Eigen::Vector2d v_Bprev = V_Bprev + w_Bprev * r_Bperp_toi;

            // The relative veolicity magnitud BEFORE collision
            const Eigen::Vector2d& n_toi = normals.row(i);

            const double vrel_prev_toi
                = (v_Aprev - v_Bprev).transpose() * n_toi;

            // solve for the impulses
            const double nr_A_toi = n_toi.transpose() * r_Aperp_toi;
            const double nr_B_toi = n_toi.transpose() * r_Bperp_toi;
            const double K = (inv_m_A + inv_m_B //
                + inv_I_A * nr_A_toi * nr_A_toi
                + inv_I_B * nr_B_toi * nr_B_toi);

            const double j
                = -1.0 / K * (1.0 + coefficient_restitution) * vrel_prev_toi;

            // update
            Eigen::Vector2d V_A_delta = inv_m_A * j * n_toi;
            Eigen::Vector2d V_B_delta = -inv_m_B * j * n_toi;
            double w_A_delta = inv_I_A * j * nr_A_toi;
            double w_B_delta = -inv_I_B * j * nr_B_toi;
            body_A.velocity.head(2) += V_A_delta;
            body_B.velocity.head(2) += V_B_delta;

            body_A.velocity(2) += w_A_delta;
            body_B.velocity(2) += w_B_delta;
        }
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

        original_ev_impacts = m_constraint_ptr->ev_impacts;
        std::sort(original_ev_impacts.begin(), original_ev_impacts.end(),
            ccd::compare_impacts_by_time<ccd::EdgeVertexImpact>);
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
