#include "distance_barrier_constraint.hpp"

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

    void DistanceBarrierConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        Eigen::MatrixXd vertices_t1 = (*vertices) + Uk;
        barriers.resize(vertices->rows() * edges->rows());
        int num_edges = int(edges->rows());
        for (int i = 0; i < int(vertices->rows()); ++i) {
            for (int j = 0; j < int(edges->rows()); ++j) {
                int vertex_a = edges->coeff(j, 0);
                int vertex_b = edges->coeff(j, 1);
                if (vertex_a == i || vertex_b == i) {
                    barriers(i * num_edges + j) = 0;
                    continue;
                }
                Eigen::VectorXd a, b, c;
                a = vertices_t1.row(vertex_a);
                b = vertices_t1.row(vertex_b);
                c = vertices_t1.row(i);

                barriers(i * num_edges + j) = distance_barrier<double>(a, b, c);
            }
        }
        /// add inf if collision
        for (auto& ev_impact : ev_impacts) {
            long i = ev_impact.vertex_index;
            long j = ev_impact.edge_index;
            barriers(i * num_edges + j)
                += std::numeric_limits<double>::infinity();
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
        for (int i = 0; i < int(vertices->rows()); ++i) {
            for (int j = 0; j < int(edges->rows()); ++j) {
                int vertex_a = edges->coeff(j, 0);
                int vertex_b = edges->coeff(j, 1);
                if (vertex_a == i || vertex_b == i) {
                    continue;
                }
                Eigen::VectorXd a, b, c;
                a = vertices_t1.row(vertex_a);
                b = vertices_t1.row(vertex_b);
                c = vertices_t1.row(i);
                Eigen::VectorXd grad = distance_barrier_grad(a, b, c);

                std::vector<int> nodes = { { vertex_a, vertex_b, i } };

                for (size_t nid = 0; nid < 3; ++nid) {
                    for (int dim = 0; dim < 2; ++dim) {
                        barriers_jacobian(
                            i * num_edges + j, nodes[nid] + num_vertices * dim)
                            = grad[2 * int(nid) + dim];
                    }
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

        Eigen::MatrixXd vertices_t1 = (*vertices) + Uk;
        barriers_hessian.clear();
        int num_vertices = int(vertices->rows());
        for (int i = 0; i < int(vertices->rows()); ++i) {
            for (int j = 0; j < int(edges->rows()); ++j) {
                Eigen::SparseMatrix<double> global_el_hessian(
                    int(vertices->size()), int(vertices->size()));

                int vertex_a = edges->coeff(j, 0);
                int vertex_b = edges->coeff(j, 1);
                if (vertex_a == i || vertex_b == i) {
                    barriers_hessian.push_back(global_el_hessian);
                    continue;
                }
                Eigen::VectorXd a, b, c;
                a = vertices_t1.row(vertex_a);
                b = vertices_t1.row(vertex_b);
                c = vertices_t1.row(i);
                Eigen::MatrixXd hess = distance_barrier_hess(a, b, c);

                std::vector<int> nodes = { { vertex_a, vertex_b, i } };
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
        }
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
        Diff::DScalar1 barrier = distance_barrier<Diff::DScalar1>(da, db, dc);

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
        Diff::DScalar2 barrier = distance_barrier<Diff::DScalar2>(da, db, dc);
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
