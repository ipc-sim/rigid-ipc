#include "rigid_body_problem.hpp"

#include <utils/flatten.hpp>
#include <utils/tensor.hpp>

#include <autodiff/finitediff.hpp>

#include <iostream>

namespace ccd {

namespace opt {

    void RigidBodyProblem2::init(const std::vector<physics::RigidBody> rbs,
        opt::BarrierConstraint& constraint)
    {
        m_constraint_ptr = std::make_shared<opt::BarrierConstraint>(constraint);
        m_assembler.init(rbs);
    }

    bool RigidBodyProblem2::simulation_step(const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        Eigen::MatrixXd x0 = m_assembler.world_vertices();
        Eigen::VectorXd g0 = m_assembler.rb_positions();

        for (auto& rb : m_assembler.m_rbs) {
            rb.position_prev = rb.position;
            rb.position += time_step * rb.velocity;
        }
        Eigen::MatrixXd x1 = m_assembler.world_vertices();

        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(x0, x1 - x0, m_assembler.m_edges,
            ev_impacts, m_constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }

    bool RigidBodyProblem2::take_step(
        const Eigen::VectorXd& rb_positions, const double time_step)
    {
        assert(m_constraint_ptr != nullptr);

        m_assembler.set_rb_positions(rb_positions);
        for (auto& rb : m_assembler.m_rbs) {
            rb.velocity = (rb.position - rb.position_prev) / time_step;
        }
        // check for collisions
        Eigen::MatrixXd x0 = m_assembler.world_vertices(/*prev=*/true);
        Eigen::MatrixXd x1 = m_assembler.world_vertices(/*prev=*/false);

        EdgeVertexImpacts ev_impacts;
        ccd::detect_edge_vertex_collisions(x0, x1 - x0, m_assembler.m_edges,
            ev_impacts, m_constraint_ptr->detection_method,
            /*reset_impacts=*/true);

        return ev_impacts.size() > 0;
    }
    void RigidBodyProblem2::update_constraint()
    {
        m_x0 = m_assembler.world_vertices(/*previous=*/true);
        m_sigma1 = m_assembler.rb_positions(/*previous=*/false);

        // base problem initial solution
        x0 = m_assembler.rb_positions(/*previous=*/true);
        num_vars = int(x0.size());
        is_dof_fixed = Eigen::MatrixXb::Zero(num_vars, 1);

        // intended displacements
        Eigen::MatrixXd x1 = m_assembler.world_vertices(/*previous=*/false);
        m_constraint_ptr->initialize(m_x0, m_assembler.m_edges, x1 - m_x0);
    }

    double RigidBodyProblem2::eval_f(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd diff = (sigma - m_sigma1);
        return 0.5 * diff.transpose() * m_assembler.m_mass_matrix * diff;
    }

    Eigen::VectorXd RigidBodyProblem2::eval_grad_f(const Eigen::VectorXd& sigma)
    {
        Eigen::VectorXd diff = (sigma - m_sigma1);
        return m_assembler.m_mass_matrix * diff;
    }

    Eigen::SparseMatrix<double> RigidBodyProblem2::eval_hessian_f_sparse(
        const Eigen::VectorXd& /*sigma*/)
    {
        return m_assembler.m_mass_matrix;
    }

    Eigen::VectorXd RigidBodyProblem2::eval_g(const Eigen::VectorXd& sigma)
    {

        Eigen::MatrixXd m_xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = m_xk - m_x0;

        if (m_constraint_ptr->update_collision_set) {
            m_constraint_ptr->detectCollisions(uk);
        }
        Eigen::VectorXd g_uk;
        m_constraint_ptr->compute_constraints(uk, g_uk);
        return g_uk;
    }

    Eigen::MatrixXd RigidBodyProblem2::eval_jac_g(const Eigen::VectorXd& sigma)
    {

        Eigen::MatrixXd m_xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = m_xk - m_x0;

        Eigen::SparseMatrix<double> jac_xk;
        m_assembler.world_vertices_gradient(sigma, jac_xk);

        if (m_constraint_ptr->update_collision_set) {
            m_constraint_ptr->detectCollisions(uk);
        }
        Eigen::MatrixXd jac_g_uk;
        m_constraint_ptr->compute_constraints_jacobian(uk, jac_g_uk);

        return jac_g_uk * jac_xk;
    }

    std::vector<Eigen::SparseMatrix<double>> RigidBodyProblem2::eval_hessian_g(
        const Eigen::VectorXd& sigma)
    {

        Eigen::MatrixXd m_xk = m_assembler.world_vertices(sigma);
        Eigen::MatrixXd uk = m_xk - m_x0;

        if (m_constraint_ptr->update_collision_set) {
            m_constraint_ptr->detectCollisions(uk);
        }

        Eigen::SparseMatrix<double> jac_xk_sigma;
        std::vector<Eigen::SparseMatrix<double>> hess_xk_sigma;
        m_assembler.world_vertices_gradient(sigma, jac_xk_sigma);
        m_assembler.world_vertices_hessian(sigma, hess_xk_sigma);

        Eigen::MatrixXd jac_g_uk;
        std::vector<Eigen::SparseMatrix<double>> hessian_g_uk;
        m_constraint_ptr->compute_constraints_jacobian(uk, jac_g_uk);
        m_constraint_ptr->compute_constraints_hessian(uk, hessian_g_uk);

        size_t num_constraints = size_t(jac_g_uk.rows());
        assert(hessian_g_uk.size() == num_constraints);
        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        gx_hessian.reserve(num_constraints);

        for (unsigned long i = 0; i < num_constraints; i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(
                jac_xk_sigma.transpose() * hessian_g_uk[i] * jac_xk_sigma
                + tensor::multiply(jac_g_uk.row(long(i)), hess_xk_sigma));
        }
        return gx_hessian;
    }

    RigidBodyProblem::RigidBodyProblem()
        : model(nullptr)
    {
    }

    RigidBodyProblem::~RigidBodyProblem() {}

    void RigidBodyProblem::initialize(
        ccd::physics::RigidBodySystem& rbs, CollisionConstraint& cstr)
    {
        model = std::make_shared<ccd ::physics::RigidBodySystem>(rbs);
        model->assemble();

        init_num_vars();
        is_dof_fixed = Eigen::MatrixXb::Zero(num_vars, 1);

        ParticlesDisplProblem::initialize(
            model->vertices, model->edges, model->displacements, cstr);
    }

    void RigidBodyProblem::init_num_vars()
    {
        num_vars = 3 * this->model->get_rigid_body_count();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Objective function and its derivatives.

    double RigidBodyProblem::eval_f(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);
        return ParticlesDisplProblem::eval_f(Uk);
    }

    Eigen::VectorXd RigidBodyProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        // TODO: compute displacement, gradient and hessian together.
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);

        Eigen::SparseMatrix<double> jac_uk_x;
        model->compute_displacements_gradient(x, jac_uk_x);
#ifdef WITH_DERIVATIVE_CHECK
        auto foo = [&](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            Eigen::MatrixXd Uk;
            this->model->compute_displacements(x, Uk);
            flatten(Uk);
            return Uk;
        };
        Eigen::MatrixXd approx_jac_uk_x;
        ccd::finite_jacobian(x, foo, approx_jac_uk_x);
        assert(compare_jacobian(Eigen::MatrixXd(jac_uk_x), approx_jac_uk_x));
#endif

        Eigen::VectorXd grad_uk = ParticlesDisplProblem::eval_grad_f(Uk);
#ifdef WITH_DERIVATIVE_CHECK
        auto bar = [&](const Eigen::VectorXd& Uk) -> double {
            return ParticlesDisplProblem::eval_f(Uk);
        };
        Eigen::VectorXd approx_grad_uk;
        ccd::finite_gradient(Uk, bar, approx_grad_uk);
        assert(compare_gradient(grad_uk, approx_grad_uk));
#endif

        return grad_uk.transpose() * jac_uk_x;
    }

    Eigen::MatrixXd RigidBodyProblem::eval_hessian_f(const Eigen::VectorXd& x)
    {
        return Eigen::MatrixXd(eval_hessian_f_sparse(x));
    }

    Eigen::SparseMatrix<double> RigidBodyProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        // TODO: compute displacement, gradient and hessian together.
        Eigen::MatrixXd Uk;
        model->compute_displacements(x, Uk);
        flatten(Uk);

        std::vector<Eigen::SparseMatrix<double>> hess_uk_x;
        model->compute_displacements_hessian(x, hess_uk_x);

        Eigen::SparseMatrix<double> jac_uk_x;
        model->compute_displacements_gradient(x, jac_uk_x);

        Eigen::SparseMatrix<double> hess_uk
            = ParticlesDisplProblem::eval_hessian_f_sparse(Uk);
        Eigen::VectorXd grad_uk = ParticlesDisplProblem::eval_grad_f(Uk);

        // chain rule
        return jac_uk_x.transpose() * hess_uk * jac_uk_x
            + tensor::multiply(grad_uk, hess_uk_x);
    }
    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Constraint function and its derivatives.

    Eigen::VectorXd RigidBodyProblem::eval_g(const Eigen::VectorXd& x)
    {
        // Compute the displacements.
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);
        // Compute the gradient of the displacements g(U(x))
        return ParticlesDisplProblem::eval_g(displacements);
    }

    // Evaluate the constraint's jacobian at point x.
    Eigen::MatrixXd RigidBodyProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        // Compute the displacements, U.
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        // Compute ∇U(x)
        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        // Compute ∇g(U)
        Eigen::MatrixXd g_gradient
            = ParticlesDisplProblem::eval_jac_g(displacements);

        // Apply the chain rule to compute ∇g(U(x)) = ∇g(U) * ∇U(x)
        return g_gradient * displacements_gradient;
    }

    // Evaluate the constraint's hessian at point x.
    std::vector<Eigen::SparseMatrix<double>> RigidBodyProblem::eval_hessian_g(
        const Eigen::VectorXd& x)
    {
        // Compute the displacements U
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        // Compute ∇U(x)
        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        // Compute ∇²U(x)
        std::vector<Eigen::SparseMatrix<double>> displacements_hessian;
        model->compute_displacements_hessian(x, displacements_hessian);

        Eigen::VectorXd _;                                  // g(U(x))
        Eigen::MatrixXd g_gradient;                         // ∇g(U)
        std::vector<Eigen::SparseMatrix<double>> g_hessian; // ∇²g(U)
        ParticlesDisplProblem::eval_g_and_gdiff(
            displacements, _, g_gradient, g_hessian);

        // Compute ∇²g(U(x))
        std::vector<Eigen::SparseMatrix<double>> gx_hessian;
        gx_hessian.reserve(g_hessian.size());
        for (unsigned long i = 0; i < g_hessian.size(); i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(displacements_gradient.transpose()
                    * g_hessian[i] * displacements_gradient
                + tensor::multiply(
                      g_gradient.row(long(i)), displacements_hessian));
        }

        return gx_hessian;
    }

    // Evaluates the constraints, jacobian, and hessian at point x.
    void RigidBodyProblem::eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& gx,
        Eigen::MatrixXd& gx_gradient,
        std::vector<Eigen::SparseMatrix<double>>& gx_hessian)
    {
        Eigen::MatrixXd displacements;
        model->compute_displacements(x, displacements);
        flatten(displacements);

        Eigen::SparseMatrix<double> displacements_gradient;
        model->compute_displacements_gradient(x, displacements_gradient);

        std::vector<Eigen::SparseMatrix<double>> displacements_hessian;
        model->compute_displacements_hessian(x, displacements_hessian);

        Eigen::MatrixXd g_gradient;
        std::vector<Eigen::SparseMatrix<double>> g_hessian;
        ParticlesDisplProblem::eval_g_and_gdiff(
            displacements, gx, g_gradient, g_hessian);

        gx_gradient = g_gradient * displacements_gradient;

        gx_hessian.clear();
        gx_hessian.reserve(g_hessian.size());
        for (unsigned long i = 0; i < g_hessian.size(); i++) {
            // Compute ∇²gᵢ(U(x)) = [∇U(x)]ᵀ * ∇²gᵢ(U) * ∇U(x) + ∇gᵢ(U) * ∇²U(x)
            gx_hessian.push_back(displacements_gradient.transpose()
                    * g_hessian[i] * displacements_gradient
                + tensor::multiply(
                      g_gradient.row(long(i)), displacements_hessian));
        }
    }
    ////////////////////////////////////////////////////////////////////////////

    // Evaluate the intermediate callback.
    bool RigidBodyProblem::eval_intermediate_callback(const Eigen::VectorXd& x)
    {
        if (intermediate_callback != nullptr) {
            Eigen::MatrixXd Uk;
            model->compute_displacements(x, Uk);
            return intermediate_callback(x, Uk);
        }

        return true;
    }

} // namespace opt
} // namespace ccd
