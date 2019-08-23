#include "distance_barrier_constraint.hpp"

#include <ccd/collision_detection.hpp>
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
        , m_barrier_epsilon(0.0)

    {
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        custom_inital_epsilon = json["custom_initial_epsilon"].get<double>();
        active_constraint_scale = json["active_constraint_scale"].get<double>();
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["custom_inital_epsilon"] = custom_inital_epsilon;
        json["active_constraint_scale"] = active_constraint_scale;

        return json;
    }

    EdgeVertexImpacts DistanceBarrierConstraint::initialize(
        const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXi& group_ids,
        const Eigen::MatrixXd& Uk)
    {
        m_barrier_epsilon = custom_inital_epsilon;
        return CollisionConstraint::initialize(vertices, edges, group_ids, Uk);
    }

    DistanceBarrierConstraint::CollisionCheck
    DistanceBarrierConstraint::get_active_barrier_set(
        const Eigen::MatrixXd& Uk, EdgeVertexCandidates& ev_barriers)
    {
        NAMED_PROFILE_POINT(
            "distance_barrier__update_collision_set", BROAD_PHASE)
        PROFILE_START(BROAD_PHASE)
        EdgeVertexCandidates ev_candidates;
        detect_edge_vertex_collision_candidates(vertices, Uk, edges, group_ids,
            ev_candidates, detection_method, m_barrier_epsilon);
        PROFILE_END(BROAD_PHASE)

        NAMED_PROFILE_POINT("distance_barrier__update_active_set", NARROW_PHASE)
        PROFILE_START(NARROW_PHASE)

        ev_barriers.clear();
        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        for (size_t i = 0; i < ev_candidates.size(); i++) {
            const EdgeVertexCandidate& ev_candidate = ev_candidates[i];
            double toi, alpha;
            bool active_impact = compute_edge_vertex_time_of_impact(
                vertices.row(edges(ev_candidate.edge_index, 0)),
                vertices.row(edges(ev_candidate.edge_index, 1)),
                vertices.row(ev_candidate.vertex_index),
                Uk.row(edges(ev_candidate.edge_index, 0)),
                Uk.row(edges(ev_candidate.edge_index, 1)),
                Uk.row(ev_candidate.vertex_index), toi, alpha);

            if (active_impact) {
                return DistanceBarrierConstraint::HAS_COLLISION;
            }

            double distance = sqrt(point_to_edge_sq_distance<double>(
                vertices_t1.row(edges(ev_candidate.edge_index, 0)),
                vertices_t1.row(edges(ev_candidate.edge_index, 1)),
                vertices_t1.row(ev_candidate.vertex_index)));

            bool distance_active
                = distance < active_constraint_scale * m_barrier_epsilon;
            if (distance_active) {
                ev_barriers.push_back(ev_candidate);
            }
        }

        return DistanceBarrierConstraint::NO_COLLISIONS;
        PROFILE_END(NARROW_PHASE)
    }

    void DistanceBarrierConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        EdgeVertexCandidates ev_candidates;
        auto check = get_active_barrier_set(Uk, ev_candidates);

        if (check == DistanceBarrierConstraint::HAS_COLLISION) {
            barriers.resize(1);
            barriers.setConstant(1, 1, std::numeric_limits<double>::infinity());
            return;
        } else {
            compute_candidates_constraints(Uk, ev_candidates, barriers);
        }
    }
    void DistanceBarrierConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian)
    {
        EdgeVertexCandidates ev_candidates;
        auto check = get_active_barrier_set(Uk, ev_candidates);
        assert(check == DistanceBarrierConstraint::NO_COLLISIONS);
        compute_candidates_constraints_jacobian(
            Uk, ev_candidates, barriers_jacobian);
    }
    void DistanceBarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        EdgeVertexCandidates ev_candidates;
        auto check = get_active_barrier_set(Uk, ev_candidates);
        assert(check == DistanceBarrierConstraint::NO_COLLISIONS);

        compute_candidates_constraints_hessian(
            Uk, ev_candidates, barriers_hessian);
    }

//    void DistanceBarrierConstraint::compute_candidates_constraints(
//        const Eigen::MatrixXd& Uk,
//        const EdgeVertexCandidates& ev_candidates,
//        Eigen::VectorXd& barriers)
//    {
//        // distance barrier is evaluated at end-positions
//        Eigen::MatrixXd vertices_t1 = vertices + Uk;

//        barriers.resize(ev_candidates.size());
//        barriers.setZero();
//        for (size_t i = 0; i < ev_candidates.size(); ++i) {
//            const auto& ev_candidate = ev_candidates[i];
//            // a and b are the endpoints of the edge; c is the vertex
//            long edge_id = ev_candidate.edge_index;
//            int a_id = edges.coeff(edge_id, 0);
//            int b_id = edges.coeff(edge_id, 1);
//            long c_id = ev_candidate.vertex_index;
//            assert(a_id != c_id && b_id != c_id);
//            Eigen::VectorXd a = vertices_t1.row(a_id);
//            Eigen::VectorXd b = vertices_t1.row(b_id);
//            Eigen::VectorXd c = vertices_t1.row(c_id);

//            barriers(int(i)) = distance_barrier<double>(a, b, c);
//        }
//    }

    void DistanceBarrierConstraint::compute_candidates_constraints_jacobian(
        const Eigen::MatrixXd& Uk,
        const EdgeVertexCandidates& ev_candidates,
        Eigen::MatrixXd& barriers_jacobian)
    {
        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        int num_vertices = int(vertices.rows());
        barriers_jacobian.resize(ev_candidates.size(), vertices.size());
        barriers_jacobian.setZero();

        for (size_t i = 0; i < ev_candidates.size(); ++i) {
            const auto& ev_candidate = ev_candidates[i];
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

    void DistanceBarrierConstraint::compute_candidates_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        const EdgeVertexCandidates& ev_candidates,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        Eigen::MatrixXd vertices_t1 = vertices + Uk;

        int num_vertices = int(vertices.rows());

        barriers_hessian.clear();
        barriers_hessian.resize(ev_candidates.size());

        for (size_t i = 0; i < ev_candidates.size(); ++i) {
            const auto& ev_candidate = ev_candidates[i];
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
                            triplets.push_back(
                                M(nodes[nid_i] + num_vertices * dim_i,
                                    nodes[nid_j] + num_vertices * dim_j,
                                    hess(2 * int(nid_i) + dim_i,
                                        2 * int(nid_j) + dim_j)));
                        }
                    }
                }
            }

            global_el_hessian.setFromTriplets(triplets.begin(), triplets.end());
            barriers_hessian[i] = global_el_hessian;
        }
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



} // namespace opt
} // namespace ccd
