#include "distance_barrier_rb_problem.hpp"

namespace ccd {

void DistanceBarrierRBProblem::update_constraints(
    const Eigen::MatrixXd& uk, const bool update_cstr_set)
{

    if (update_cstr_set) {
        constraint_.update_collision_set(uk);
    }
    constraint_.update_active_set(uk);
}

Eigen::VectorXd DistanceBarrierRBProblem::eval_g(const Eigen::VectorXd& sigma)
{
    return eval_g_(sigma, true);
}

Eigen::VectorXd DistanceBarrierRBProblem::eval_g_(
    const Eigen::VectorXd& sigma, const bool update_cstr_set)
{

    Eigen::MatrixXd xk = m_assembler.world_vertices(sigma);
    Eigen::MatrixXd uk = xk - vertices_t0;

    update_constraints(xk, update_cstr_set);

    Eigen::VectorXd g_uk;
    constraint_.compute_constraints(uk, g_uk);
    return g_uk;
}

/// @brief eval_jac_g evaluates constraints jacobian at point x
Eigen::MatrixXd DistanceBarrierRBProblem::eval_jac_g(
    const Eigen::VectorXd& sigma)
{

    Eigen::MatrixXd xk = m_assembler.world_vertices(sigma);
    Eigen::MatrixXd uk = xk - vertices_t0;
    update_constraints(xk, true);

    Eigen::MatrixXd jac_g;
    jac_g.resize(constraint_.m_num_active_constraints, num_vars_);
    jac_g.setZero();

    typedef ccd::DistanceBarrierDiff Diff;
    Diff::activate();
    Diff::D1Vector3d position_E;
    Diff::D1Vector3d position_V;
    for (size_t i = 0; i < constraint_.m_ev_distance_active.size(); ++i) {
        const auto& ev_candidate = constraint_.m_ev_distance_active[i];

        // get differentiable variables for the position of the rigid bodies
        const long v_id = ev_candidate.vertex_index;
        const int e0_id = m_assembler.m_edges.coeff(ev_candidate.edge_index, 0);
        const int e1_id = m_assembler.m_edges.coeff(ev_candidate.edge_index, 1);
        int body_V_id, lv_id, body_E_id, le0_id, le1_id;

        m_assembler.global_to_local(v_id, body_V_id, lv_id);
        m_assembler.global_to_local(e0_id, body_E_id, le0_id);
        m_assembler.global_to_local(e1_id, body_E_id, le1_id);

        position_V = ccd::DistanceBarrierDiff::d1vars(
            0, sigma.segment(3 * body_V_id, 3));
        position_E = ccd::DistanceBarrierDiff::d1vars(
            3, sigma.segment(3 * body_E_id, 3));

        // get position of the vertices involved in the collision wrt the
        // variables
        Diff::D1Vector3d da = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_E, le0_id);
        Diff::D1Vector3d db = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_E, le1_id);
        Diff::D1Vector3d dc = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_V, lv_id);

        Diff::DDouble1 barrier
            = constraint_.distance_barrier<Diff::DDouble1>(da, db, dc);
        Eigen::VectorXd gradient = barrier.getGradient();

        jac_g.row(int(i)).segment(3 * body_V_id, 3) = gradient.segment(0, 3);
        jac_g.row(int(i)).segment(3 * body_E_id, 3) = gradient.segment(3, 3);
    }
    return jac_g;
}

std::vector<Eigen::SparseMatrix<double>>
DistanceBarrierRBProblem::eval_hessian_g(const Eigen::VectorXd& sigma)
{
    std::vector<Eigen::SparseMatrix<double>> gx_hessian;
    gx_hessian.reserve(size_t(constraint_.m_num_active_constraints));

    typedef ccd::DistanceBarrierDiff Diff;
    Diff::activate();
    Diff::D2Vector3d position_E;
    Diff::D2Vector3d position_V;

    typedef Eigen::Triplet<double> M;
    std::vector<M> triplets;

    for (size_t i = 0; i < constraint_.m_ev_distance_active.size(); ++i) {
        const auto& ev_candidate = constraint_.m_ev_distance_active[i];
        // get differentiable variables for the position of the rigid bodies
        const long v_id = ev_candidate.vertex_index;
        const int e0_id = m_assembler.m_edges.coeff(ev_candidate.edge_index, 0);
        const int e1_id = m_assembler.m_edges.coeff(ev_candidate.edge_index, 1);
        int body_V_id, lv_id, body_E_id, le0_id, le1_id;

        m_assembler.global_to_local(v_id, body_V_id, lv_id);
        m_assembler.global_to_local(e0_id, body_E_id, le0_id);
        m_assembler.global_to_local(e1_id, body_E_id, le1_id);

        position_V = ccd::DistanceBarrierDiff::d2vars(
            0, sigma.segment(3 * body_V_id, 3));
        position_E = ccd::DistanceBarrierDiff::d2vars(
            3, sigma.segment(3 * body_E_id, 3));

        // get position of the vertices involved in the collision wrt the
        // variables
        Diff::D2Vector3d da = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_E, le0_id);
        Diff::D2Vector3d db = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_E, le1_id);
        Diff::D2Vector3d dc = m_assembler.m_rbs[size_t(body_E_id)].world_vertex(
            position_V, lv_id);

        Diff::DDouble2 barrier
            = constraint_.distance_barrier<Diff::DDouble2>(da, db, dc);
        Eigen::MatrixXd hessian = barrier.getHessian();

        triplets.clear();
        triplets.reserve(6 * 6);
        int bodies[2] = { body_V_id, body_E_id };

        for (int b_i = 0; b_i < 2; b_i++) {
            for (int b_j = 0; b_j < 2; b_j++) {
                for (int dim_i = 0; dim_i < 3; dim_i++) {
                    for (int dim_j = 0; dim_j < 3; dim_j++) {
                        double v = hessian(3 * b_i + dim_i, 3 * b_j + dim_j);
                        int r = 3 * bodies[b_i] + dim_i;
                        int c = 3 * bodies[b_j] + dim_j;
                        triplets.push_back(M(r, c, v));
                    }
                }
            }
        }
    }
    /// fill in for collisions
    for (size_t i = gx_hessian.size();
         i < size_t(constraint_.m_num_active_constraints); ++i) {
        Eigen::SparseMatrix<double> global_el_hessian(num_vars_, num_vars_);
        gx_hessian.push_back(global_el_hessian);
    }

    return gx_hessian;
}

} // namespace ccd
