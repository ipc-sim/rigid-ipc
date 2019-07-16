#include "collision_constraint.hpp"

#include <ccd/prune_impacts.hpp>
#include <ccd/time_of_impact.hpp>

#include <profiler.hpp>

namespace ccd {
namespace opt {

    CollisionConstraint::CollisionConstraint(const std::string& name)
        : detection_method(HASH_GRID)
        , extend_collision_set(true)
        , name_(name)
    {
    }

    CollisionConstraint::~CollisionConstraint() {}

    void CollisionConstraint::initialize(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2i& E,
        const Eigen::MatrixXd& Uk)
    {
        vertices = std::make_shared<const Eigen::MatrixX2d>(V);
        edges = std::make_shared<const Eigen::MatrixX2i>(E);

        ev_impacts.clear();
        ee_impacts.clear();
        edge_impact_map.resize(E.rows());
        edge_impact_map.setZero();
        detectCollisions(Uk);
    }

    void CollisionConstraint::detectCollisions(const Eigen::MatrixXd& Uk)
    {
        edge_impact_map.resize(edges->rows());
        ccd::detect_edge_vertex_collisions(*vertices, Uk, *edges, ev_impacts,
            detection_method, /*reset_impacts=*/!extend_collision_set);
        ccd::convert_edge_vertex_to_edge_edge_impacts(
            *edges, ev_impacts, ee_impacts);
        num_pruned_impacts = prune_impacts(ee_impacts, edge_impact_map);
    }

    // -------------------------------------------------------------------
    // Assembly of global Matrices
    // -------------------------------------------------------------------

    void CollisionConstraint::assemble_hessian(
        const std::vector<DScalar>& constraints,
        std::vector<Eigen::SparseMatrix<double>>& hessian)
    {
        int num_vertices = int(vertices->rows());

        hessian.clear();
        hessian.reserve(ee_impacts.size() * 2);

        std::vector<Eigen::Triplet<double>> coefficients;
        coefficients.reserve(64);

        for (size_t ee = 0; ee < ee_impacts.size(); ++ee) {
            int nodes[4];
            get_impact_nodes(ee_impacts[ee], nodes);

            // Note: global gradient is sorted as x,x,x,...y,y,y
            // while local gradient is sorted as x,y,x,y,...,x,y

            for (size_t k = 0; k < 2; k++) {
                Matrix8d local_hessian = constraints[2 * ee + k].getHessian();

                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        coefficients.push_back(Eigen::Triplet<double>(
                            nodes[i], nodes[j], local_hessian(2 * i, 2 * j)));
                        coefficients.push_back(
                            Eigen::Triplet<double>(nodes[i] + num_vertices,
                                nodes[j], local_hessian(2 * i + 1, 2 * j)));
                        coefficients.push_back(Eigen::Triplet<double>(
                            nodes[i] + num_vertices, nodes[j] + num_vertices,
                            local_hessian(2 * i + 1, 2 * j + 1)));
                        coefficients.push_back(Eigen::Triplet<double>(nodes[i],
                            nodes[j] + num_vertices,
                            local_hessian(2 * i, 2 * j + 1)));
                    }
                }

                Eigen::SparseMatrix<double> global_hessian(
                    int(vertices->size()), int(vertices->size()));
                global_hessian.setFromTriplets(
                    coefficients.begin(), coefficients.end());
                hessian.push_back(global_hessian);
                coefficients.clear();
            }
        }
    }

    void CollisionConstraint::assemble_jacobian(
        const std::vector<DScalar>& constraints, Eigen::MatrixXd& jacobian)
    {
        std::vector<DoubleTriplet> tripletList;
        assemble_jacobian_triplets(constraints, tripletList);

        jacobian.resize(int(constraints.size()), vertices->size());
        jacobian.setZero();

        for (auto& triplet : tripletList) {
            jacobian(triplet.row(), triplet.col()) = triplet.value();
        }
    }

    void CollisionConstraint::assemble_jacobian(
        const std::vector<DScalar>& constraints,
        Eigen::SparseMatrix<double>& jacobian)
    {
        std::vector<DoubleTriplet> tripletList;
        assemble_jacobian_triplets(constraints, tripletList);
        jacobian.resize(int(constraints.size()), int(vertices->size()));
        jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void CollisionConstraint::assemble_jacobian_triplets(
        const std::vector<DScalar>& constraints,
        std::vector<DoubleTriplet>& tripletList)
    {
        tripletList.clear();
        tripletList.reserve(2 * ee_impacts.size());

        const int num_vertices = int(vertices->rows());

        for (size_t ee = 0; ee < ee_impacts.size(); ++ee) {
            int nodes[4];
            get_impact_nodes(ee_impacts[ee], nodes);

            for (size_t k = 0; k < 2; k++) {
                Vector8d local_jac = constraints[2 * ee + k].getGradient();

                for (int i = 0; i < 4; i++) {
                    tripletList.push_back(DoubleTriplet(
                        int(2 * ee + k), nodes[i], local_jac(2 * i + 0)));
                    tripletList.push_back(DoubleTriplet(int(2 * ee + k),
                        nodes[i] + num_vertices, local_jac(2 * i + 1)));
                }
            }
        }
    }

    void CollisionConstraint::assemble_constraints(
        const std::vector<DScalar>& constraints, Eigen::VectorXd& g_uk)
    {
        g_uk.resize(int(constraints.size()));
        for (size_t i = 0; i < constraints.size(); ++i) {
            g_uk(int(i)) = constraints[i].getValue();
        }
    }
    // Helper Functions
    //  -------------------------------------

    void CollisionConstraint::get_impact_nodes(
        const EdgeEdgeImpact& ee_impact, int nodes[4])
    {
        Eigen::Vector2i e_ij = edges->row(ee_impact.impacted_edge_index);
        Eigen::Vector2i e_kl = edges->row(ee_impact.impacting_edge_index);

        nodes[0] = e_ij(0);
        nodes[1] = e_ij(1);
        nodes[2] = e_kl(0);
        nodes[3] = e_kl(1);
    }

    template <typename T>
    bool CollisionConstraint::compute_toi_alpha(
        const ImpactTData<T>& data, T& toi, T& alpha_ij, T& alpha_kl)
    {

        bool success;
        toi = T(0);
        alpha_ij = T(0);
        alpha_kl = T(data.impacting_side);

        // clang-format off
        success = ccd::autodiff::compute_edge_vertex_time_of_impact<T>(
                    data.v[0], data.v[1], data.v[data.impacting_index()],
                    data.u[0], data.u[1], data.u[data.impacting_index()], toi);
        success = success && ccd::autodiff::temporal_parameterization_to_spatial<T>(
                    data.v[0], data.v[1], data.v[data.impacting_index()],
                    data.u[0], data.u[1], data.u[data.impacting_index()], toi, alpha_ij);
        // clang-format on

        return success;
    }

    template <>
    ImpactTData<double> CollisionConstraint::get_impact_data<double>(
        const Eigen::MatrixXd& displacements, const EdgeEdgeImpact ee_impact)
    {
        Eigen::Vector2i e_ij = edges->row(ee_impact.impacted_edge_index);
        Eigen::Vector2i e_kl = edges->row(ee_impact.impacting_edge_index);

        ImpactTData<double> data;
        slice_vector(*vertices, e_ij, e_kl, data.v);
        slice_vector(displacements, e_ij, e_kl, data.u);
        data.impacting_side = ee_impact.impacting_node();
        return data;
    }

    template <>
    ImpactTData<DScalar> CollisionConstraint::get_impact_data<DScalar>(
        const Eigen::MatrixXd& displacements, const EdgeEdgeImpact ee_impact)
    {
        ImpactTData<double> data
            = get_impact_data<double>(displacements, ee_impact);

        ImpactTData<DScalar> diff_data;
        diff_data.v = data.v;
        diff_data.impacting_side = data.impacting_side;

        diff_data.u[0] = dvector(0, data.u[0]);
        diff_data.u[1] = dvector(2, data.u[1]);
        diff_data.u[2] = dvector(4, data.u[2]);
        diff_data.u[3] = dvector(6, data.u[3]);
        return diff_data;
    }

    void slice_vector(const Eigen::MatrixX2d& data,
        const Eigen::Vector2i e_ij,
        Eigen::Vector2i e_kl,
        std::array<Eigen::Vector2d, 4>& d)
    {
        d[0] = data.row(e_ij(0));
        d[1] = data.row(e_ij(1));
        d[2] = data.row(e_kl(0));
        d[3] = data.row(e_kl(1));
    }

    template bool CollisionConstraint::compute_toi_alpha<double>(
        const ImpactTData<double>& data,
        double& toi,
        double& alpha_ij,
        double& aplha_kl);

    template bool CollisionConstraint::compute_toi_alpha<DScalar>(
        const ImpactTData<DScalar>& data,
        DScalar& toi,
        DScalar& alpha_ij,
        DScalar& aplha_kl);

    // Check if a type, T, is differentiable (differentiable<T>())
    template <> bool differentiable<DScalar>() { return true; }
    template <> bool differentiable<double>() { return false; }
} // namespace opt
} // namespace ccd
