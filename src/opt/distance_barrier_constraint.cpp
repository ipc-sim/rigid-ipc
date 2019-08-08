#include "distance_barrier_constraint.hpp"

#include <igl/slice_mask.h>

#include <opt/barrier.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    DistanceBarrierConstraint::DistanceBarrierConstraint()
        : DistanceBarrierConstraint("distance_barrier_constraint")
    {
    }

    DistanceBarrierConstraint::DistanceBarrierConstraint(
        const std::string& name)
        : CollisionConstraint(name)
        , custom_inital_epsilon(1.0)
        , active_constraint_scale(1.5)
        , use_distance_hashgrid(true)
        , m_barrier_epsilon(0.0)
        , m_num_active_constraints(0)

    {
        extend_collision_set = false;
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        custom_inital_epsilon = json["custom_initial_epsilon"].get<double>();
        use_distance_hashgrid = json["use_distance_hashgrid"].get<bool>();
        active_constraint_scale = json["active_constraint_scale"].get<double>();
        extend_collision_set = false;
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["custom_inital_epsilon"] = custom_inital_epsilon;
        json["use_distance_hashgrid"] = use_distance_hashgrid;
        json["active_constraint_scale"] = active_constraint_scale;

        return json;
    }

    void DistanceBarrierConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXi& group_ids,
        const Eigen::MatrixXd& Uk)
    {
        m_barrier_epsilon = custom_inital_epsilon;
        CollisionConstraint::initialize(vertices, edges, group_ids, Uk);
    }

    int DistanceBarrierConstraint::number_of_constraints()
    {
        return m_num_active_constraints;
    }

    void DistanceBarrierConstraint::detectCollisions(const Eigen::MatrixXd& Uk)
    {
        // We want the distance barrier and the CCD to have a common broad-phase
        EdgeVertexCandidates ev_candidates;

        PROFILE_POINT("distance_evaluation")
        NAMED_PROFILE_POINT("distance_evaluation__broad_phase", BROAD_PHASE)
        NAMED_PROFILE_POINT("distance_evaluation__narrow_phase", NARROW_PHASE)

        PROFILE_START()
        PROFILE_START(BROAD_PHASE)
        detect_edge_vertex_collision_candidates(vertices, Uk, edges, group_ids,
            ev_candidates, detection_method, m_barrier_epsilon);
        PROFILE_END(BROAD_PHASE)

        PROFILE_START(NARROW_PHASE)
        // Active candidates are those that are within the epsilon bound.
        Eigen::MatrixXd vertices_t1 = vertices + Uk;
        filter_active_barriers(vertices_t1, ev_candidates,
            active_constraint_scale, m_ev_distance_active);

        m_constraint_map.resize(vertices.rows() * edges.rows());
        m_constraint_map.setConstant(-1);
        m_num_active_constraints = int(m_ev_distance_active.size());

        for (size_t i = 0; i < m_ev_distance_active.size(); ++i) {
            const auto& ev_candidate = ev_candidates[i];
            long edge_id = ev_candidate.edge_index;
            long vertex_id = ev_candidate.vertex_index;
            m_constraint_map(vertex_id * edges.rows() + edge_id) = int(i);
        }

        // Do the CCD narrow-phase to get the actual impacts.
        ccd::detect_edge_vertex_collisions_from_candidates(vertices, Uk, edges,
            ev_candidates, ev_impacts, /*reset_impacts=*/!extend_collision_set);

        // Extend to include collisions.
        for (auto& ev_impact : ev_impacts) {
            long vertex_id = ev_impact.vertex_index;
            long edge_id = ev_impact.edge_index;
            long c_id = vertex_id * edges.rows() + edge_id;

            if (m_constraint_map(c_id) == -1) {
                m_constraint_map(c_id) = m_num_active_constraints;
                m_num_active_constraints += 1;
            }
        }
        PROFILE_END(NARROW_PHASE)
        PROFILE_END()
    }

    void DistanceBarrierConstraint::filter_active_barriers(
        const Eigen::MatrixXd& vertices,
        const EdgeVertexCandidates& candidates,
        const double thr,
        EdgeVertexCandidates& active)
    {
        active.clear();
        for (size_t i = 0; i < candidates.size(); ++i) {
            const auto& ev_candidate = candidates[i];
            // a and b are the endpoints of the edge; c is the vertex
            long edge_id = ev_candidate.edge_index;
            int a_id = edges.coeff(edge_id, 0);
            int b_id = edges.coeff(edge_id, 1);
            long c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices.row(a_id);
            Eigen::VectorXd b = vertices.row(b_id);
            Eigen::VectorXd c = vertices.row(c_id);

            double distance = sqrt(point_to_edge_sq_distance<double>(a, b, c));
            bool is_active = distance < thr * m_barrier_epsilon;
            if (is_active) {
                active.push_back(ev_candidate);
            }
        }
    }

    void DistanceBarrierConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        // distance barrier is evaluated at end-positions
        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        barriers.resize(m_num_active_constraints);
        barriers.setZero();
        for (size_t i = 0; i < m_ev_distance_active.size(); ++i) {
            const auto& ev_candidate = m_ev_distance_active[i];
            // a and b are the endpoints of the edge; c is the vertex
            long edge_id = ev_candidate.edge_index;
            int a_id = edges.coeff(edge_id, 0);
            int b_id = edges.coeff(edge_id, 1);
            long c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices_t1.row(a_id);
            Eigen::VectorXd b = vertices_t1.row(b_id);
            Eigen::VectorXd c = vertices_t1.row(c_id);

            barriers(int(i)) = distance_barrier<double>(a, b, c);
        }

        /// add inf if collision
        for (auto& ev_impact : ev_impacts) {
            long vertex_id = ev_impact.vertex_index;
            long edge_id = ev_impact.edge_index;
            int ac_id = m_constraint_map(vertex_id * edges.rows() + edge_id);
            barriers(ac_id) += std::numeric_limits<double>::infinity();
        }
    }

    void DistanceBarrierConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian)
    {
        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        int num_vertices = int(vertices.rows());
        barriers_jacobian.resize(m_num_active_constraints, vertices.size());
        barriers_jacobian.setZero();

        for (size_t i = 0; i < m_ev_distance_active.size(); ++i) {
            const auto& ev_candidate = m_ev_distance_active[i];
            long edge_id = ev_candidate.edge_index;
            int a_id = edges.coeff(edge_id, 0);
            int b_id = edges.coeff(edge_id, 1);
            int c_id = int(ev_candidate.vertex_index);
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices_t1.row(a_id);
            Eigen::VectorXd b = vertices_t1.row(b_id);
            Eigen::VectorXd c = vertices_t1.row(c_id);

            Eigen::VectorXd grad = distance_barrier_grad(a, b, c);

            std::vector<int> nodes = { { a_id, b_id, c_id } };

            for (size_t nid = 0; nid < 3; ++nid) {
                for (int dim = 0; dim < 2; ++dim) {
                    barriers_jacobian(int(i), nodes[nid] + num_vertices * dim)
                        = grad[2 * int(nid) + dim];
                }
            }
        }
    }

    void DistanceBarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        int num_vertices = int(vertices.rows());

        barriers_hessian.clear();
        barriers_hessian.reserve(size_t(m_num_active_constraints));
        for (size_t i = 0; i < m_ev_distance_active.size(); ++i) {
            const auto& ev_candidate = m_ev_distance_active[i];
            Eigen::SparseMatrix<double> global_el_hessian(
                int(vertices.size()), int(vertices.size()));

            long edge_id = ev_candidate.edge_index;
            int a_id = edges.coeff(edge_id, 0);
            int b_id = edges.coeff(edge_id, 1);
            int c_id = int(ev_candidate.vertex_index);
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices_t1.row(a_id);
            Eigen::VectorXd b = vertices_t1.row(b_id);
            Eigen::VectorXd c = vertices_t1.row(c_id);

            Eigen::MatrixXd hess = distance_barrier_hess(a, b, c);

            std::vector<int> nodes = { { a_id, b_id, c_id } };
            triplets.clear();
            triplets.reserve(3 * 3 * 2 * 2);
            for (size_t nid_i = 0; nid_i < 3; ++nid_i) {
                for (size_t nid_j = 0; nid_j < 3; ++nid_j) {
                    for (int dim_i = 0; dim_i < 2; ++dim_i) {
                        for (int dim_j = 0; dim_j < 2; ++dim_j) {
                            M(nodes[nid_i] + num_vertices * dim_i,
                                nodes[nid_j] + num_vertices * dim_j,
                                hess(2 * int(nid_i) + dim_i,
                                    2 * int(nid_j) + dim_j));
                        }
                    }
                }
            }

            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());
            barriers_hessian.push_back(global_el_hessian);
        }

        /// fill in for collisions
        for (size_t i = barriers_hessian.size();
             i < size_t(m_num_active_constraints); ++i) {
            Eigen::SparseMatrix<double> global_el_hessian(
                int(vertices.size()), int(vertices.size()));
            barriers_hessian.push_back(global_el_hessian);
        }
        assert(int(barriers_hessian.size()) == m_num_active_constraints);
    }

    void DistanceBarrierConstraint::compute_constraints_and_derivatives(
        const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& barriers,
        Eigen::MatrixXd& barriers_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        compute_constraints(Uk, barriers);
        compute_constraints_jacobian(Uk, barriers_jacobian);
        compute_constraints_hessian(Uk, barriers_hessian);
    }

    template <typename T>
    T DistanceBarrierConstraint::distance_barrier(
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c)
    {
        T distance = sqrt(point_to_edge_sq_distance<T>(a, b, c));
        return opt::spline_barrier<T>(distance, m_barrier_epsilon);
    }

    Eigen::VectorXd DistanceBarrierConstraint::distance_barrier_grad(
        const Eigen::VectorXd& a,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& c)
    {
        Diff::activate();
        Diff::D1VectorXd da = Diff::d1vars(0, a);
        Diff::D1VectorXd db = Diff::d1vars(2, b);
        Diff::D1VectorXd dc = Diff::d1vars(4, c);
        Diff::DDouble1 barrier = distance_barrier<Diff::DDouble1>(da, db, dc);

        return barrier.getGradient();
    }

    Eigen::MatrixXd DistanceBarrierConstraint::distance_barrier_hess(
        const Eigen::VectorXd& a,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& c)
    {
        Diff::activate();
        Diff::D2VectorXd da = Diff::d2vars(0, a);
        Diff::D2VectorXd db = Diff::d2vars(2, b);
        Diff::D2VectorXd dc = Diff::d2vars(4, c);
        Diff::DDouble2 barrier = distance_barrier<Diff::DDouble2>(da, db, dc);
        return barrier.getHessian();
    }

    template <typename T>
    T point_to_edge_sq_distance(const Eigen::Matrix<T, Eigen::Dynamic, 1>& a,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& b,
        const Eigen::Matrix<T, Eigen::Dynamic, 1>& c)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> ab = b - a;
        Eigen::Matrix<T, Eigen::Dynamic, 1> ac = c - a;
        Eigen::Matrix<T, Eigen::Dynamic, 1> bc = c - b;

        T e = ac.dot(ab);

        // Handle cases where c projects outside ab
        if (e <= T(0.0)) {
            return ac.dot(ac);
        }
        T f = ab.dot(ab);
        if (e >= f) {
            return bc.dot(bc);
        }

        // Handle cases where c projects onto ab
        return ac.dot(ac) - e * e / f;
    }

} // namespace opt
} // namespace ccd
