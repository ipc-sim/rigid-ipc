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
        , barrier_epsilon(0.0)
    {
        extend_collision_set = false;
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        custom_inital_epsilon = json["custom_initial_epsilon"].get<double>();
        extend_collision_set = false;
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["custom_inital_epsilon"] = custom_inital_epsilon;

        return json;
    }

    void DistanceBarrierConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk)
    {
        this->barrier_epsilon = 0.0; // Temporary until collisions are computed.
        CollisionConstraint::initialize(vertices, edges, Uk);
        this->barrier_epsilon = custom_inital_epsilon;
    }

    int DistanceBarrierConstraint::number_of_constraints()
    {
        return int(vertices->rows() * edges->rows());
    }

    void DistanceBarrierConstraint::compute_candidate_intersections_brute_force(
        const Eigen::MatrixXd& vertices_t1,
        EdgeVertexCandidates& ev_candidates){
        detect_edge_vertex_collisions_brute_force(
            vertices_t1, Eigen::MatrixXd::Zero(vertices_t1.rows(),
            vertices_t1.cols()), *edges, ev_candidates);
    }

    void DistanceBarrierConstraint::compute_candidate_intersections_hashgrid(
        const Eigen::MatrixXd& vertices_t1,
        EdgeVertexCandidates& ev_candidates){
        ev_candidates.clear();

        HashGrid hashgrid;
        Eigen::MatrixXd displacements = Eigen::MatrixXd::Zero(
            vertices_t1.rows(), vertices_t1.cols());
        hashgrid.resize(vertices_t1, displacements, *edges, barrier_epsilon);
        hashgrid.addVertices(vertices_t1, displacements, barrier_epsilon);
        hashgrid.addEdges(vertices_t1, displacements, *edges, barrier_epsilon);

        // Assume checking if vertex is and end-point of the edge is done by
        // `hashgrid.getVertexEdgePairs(...)`.
        hashgrid.getVertexEdgePairs(*edges, ev_candidates);
    }

    void DistanceBarrierConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        Eigen::MatrixXd vertices_t1 = (*vertices) + Uk;
        barriers = Eigen::VectorXd::Zero(vertices->rows() * edges->rows());

        int num_edges = int(edges->rows());

        EdgeVertexCandidates ev_candidates;
        compute_candidate_intersections_hashgrid(vertices_t1, ev_candidates);
        for(const auto& ev_candidate : ev_candidates){
            // a and b are the endpoints of the edge; c is the vertex
            int edge_id = ev_candidate.edge_index;
            int a_id = edges->coeff(edge_id, 0);
            int b_id = edges->coeff(edge_id, 1);
            int c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices_t1.row(a_id);
            Eigen::VectorXd b = vertices_t1.row(b_id);
            Eigen::VectorXd c = vertices_t1.row(c_id);

            barriers(c_id * num_edges + edge_id) =
                distance_barrier<double>(a, b, c);
        }

        /// add inf if collision
        for (auto& ev_impact : ev_impacts) {
            long i = ev_impact.vertex_index;
            long j = ev_impact.edge_index;
            barriers(i * num_edges + j)
                += std::numeric_limits<double>::infinity();
        }
    }

    void DistanceBarrierConstraint::compute_active_mask(
        const Eigen::MatrixXd& Uk, Eigen::VectorXb& mask)
    {
        Eigen::VectorXd barriers;
        compute_constraints(Uk, barriers);
        mask.resize(barriers.size());
        for (long i = 0; i < mask.size(); i++) {
            mask(i) = barriers(i) != 0.0;
        }
    }

    void DistanceBarrierConstraint::compute_constraints_jacobian(
        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian)
    {
        Eigen::MatrixXd vertices_t1 = (*vertices) + Uk;
        barriers_jacobian.resize(
            vertices->rows() * edges->rows(), vertices->size());
        barriers_jacobian.setZero();

        int num_vertices = int(vertices->rows());
        int num_edges = int(edges->rows());

        EdgeVertexCandidates ev_candidates;
        compute_candidate_intersections_hashgrid(vertices_t1, ev_candidates);
        for(const auto& ev_candidate : ev_candidates){
            int edge_id = ev_candidate.edge_index;
            int a_id = edges->coeff(edge_id, 0);
            int b_id = edges->coeff(edge_id, 1);
            int c_id = ev_candidate.vertex_index;
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorXd a = vertices_t1.row(a_id);
            Eigen::VectorXd b = vertices_t1.row(b_id);
            Eigen::VectorXd c = vertices_t1.row(c_id);

            Eigen::VectorXd grad = distance_barrier_grad(a, b, c);

            std::vector<int> nodes = { { a_id, b_id, c_id } };

            for (size_t nid = 0; nid < 3; ++nid) {
                for (int dim = 0; dim < 2; ++dim) {
                    barriers_jacobian(
                        c_id * num_edges + edge_id,
                        nodes[nid] + num_vertices * dim)
                        = grad[2 * int(nid) + dim];
                }
            }
        }

        // TODO: This code removes unactive dof, but it breaks the tests
        Eigen::VectorXb mask;
        compute_active_mask(Uk, mask);
        Eigen::MatrixXd active_barriers_jacobian;
        assert(mask.size() == barriers_jacobian.rows());
        igl::slice_mask(barriers_jacobian, mask,
            Eigen::VectorXb::Ones(Uk.size()), active_barriers_jacobian);
        barriers_jacobian = active_barriers_jacobian;
    }

    void DistanceBarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        typedef Eigen::Triplet<double> M;
        std::vector<M> triplets;

        Eigen::MatrixXd vertices_t1 = (*vertices) + Uk;
        barriers_hessian.clear();
        int num_vertices = int(vertices->rows());

        EdgeVertexCandidates ev_candidates;
        compute_candidate_intersections_hashgrid(vertices_t1, ev_candidates);
        for(const auto& ev_candidate : ev_candidates){
            Eigen::SparseMatrix<double> global_el_hessian(
                int(vertices->size()), int(vertices->size()));

            int edge_id = ev_candidate.edge_index;
            int a_id = edges->coeff(edge_id, 0);
            int b_id = edges->coeff(edge_id, 1);
            int c_id = ev_candidate.vertex_index;
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
                                hess(2 * nid_i + dim_i, 2 * nid_j + dim_j));
                        }
                    }
                }
            }

            global_el_hessian.setFromTriplets(
                triplets.begin(), triplets.end());
                barriers_hessian.push_back(global_el_hessian);
        }

        // TODO: This code removes unactive dof, but it breaks the tests
        Eigen::VectorXb mask;
        compute_active_mask(Uk, mask);
        std::vector<Eigen::SparseMatrix<double>> active_barriers_hessian;
        assert(mask.size() == barriers_hessian.size());
        for (int i = 0; i < mask.size(); i++) {
            if (mask(i)) {
                active_barriers_hessian.push_back(barriers_hessian[i]);
            }
        }
        barriers_hessian = active_barriers_hessian;
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
        return opt::spline_barrier<T>(distance, barrier_epsilon);
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
