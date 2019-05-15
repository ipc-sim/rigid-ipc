#include <ccd/collision_volume_diff.hpp>

#include <autodiff/finitediff.hpp>
#include <autogen/collision_volume.hpp>
#include <ccd/time_of_impact.hpp>

#include <iostream>

namespace ccd {
namespace autodiff {

    template <typename T>
    T collision_volume(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const ImpactNode impact_node,
        const double epsilon)
    {
        T toi, alpha;
        bool success;

        switch (impact_node) {
        case vK:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vk, Ui, Uj, Uk, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vk, Ui, Uj, Uk, toi, alpha);
            break;
        case vL:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vl, Ui, Uj, Ul, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vl, Ui, Uj, Ul, toi, alpha);
            break;
        case vI:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vi, Uk, Ul, Ui, toi);
            alpha = T(0);
            break;
        case vJ:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vj, Uk, Ul, Uj, toi);
            alpha = T(1);
            break;
        }

        if (!success) {
            return T(0.0);
        }

        return ccd::autogen::space_time_collision_volume(
            Vi, Vj, Ui, Uj, toi, alpha, epsilon);
    }

    // -----------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Volumes Derivatives
    // -----------------------------------------------------------------------------

    void compute_volumes_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::VectorXd& volumes)
    {
        compute_constraints_dense_refresh_toi(V, U, E, ee_impacts, edge_impact_map,
            epsilon, collision_volume<double>, volumes);
    }

    void compute_volumes_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        Eigen::MatrixXd& volume_grad)
    {
        compute_constraints_dense_gradient(V, U, E, ee_impacts, edge_impact_map,
            epsilon, collision_volume<DScalar>, volume_grad);
    }

    void compute_volumes_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        std::vector<Eigen::SparseMatrix<double>>& volume_hessian)
    {
        // TODO: this should be dense too
        compute_constraints_hessian(V, U, E, ee_impacts, edge_impact_map,
            epsilon, collision_volume<DScalar>, volume_hessian);
    }

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Volumes & Derivatives
    // -----------------------------------------------------------------------------
    double collision_volume_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon)
    {
        return collision_constraint_refresh_toi(vertices, displacements, edges,
            impact, edge_id, epsilon, collision_volume<double>);
    }

    void collision_volume_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        Eigen::VectorXd& gradient)
    {
        Eigen::SparseMatrix<double> sparse_gradient = gradient.sparseView();
        collision_constraint_grad(vertices, displacements, edges, impact,
            edge_id, epsilon, collision_volume<DScalar>, sparse_gradient);
        gradient = Eigen::VectorXd(sparse_gradient);
    }

    void collision_volume_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        Eigen::MatrixXd& hessian)
    {
        Eigen::SparseMatrix<double> sparse_hessian = hessian.sparseView();
        collision_constraint_hessian(vertices, displacements, edges, impact,
            edge_id, epsilon, collision_volume<DScalar>, sparse_hessian);
        hessian = Eigen::MatrixXd(sparse_hessian);
    }

    template <int I>
    void collision_volume_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const int order, Eigen::Matrix<double, Eigen::Dynamic, I>& derivative)
    {
        Eigen::SparseMatrix<double> sparse_derivative = derivative.sparseView();
        collision_constraint_derivative(vertices, displacements, edges, impact,
            edge_id, epsilon, order, collision_volume<DScalar>,
            sparse_derivative);
        derivative
            = Eigen::Matrix<double, Eigen::Dynamic, I>(sparse_derivative);
    }

    //-------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Volumes Derivatives
    //-------------------------------------------------------------------------
    DScalar collision_volume_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon)
    {
        return collision_constraint_differentiable(Vi, Vj, Vk, Vl, Ui, Uj, Uk,
            Ul, impact_node, epsilon, collision_volume<DScalar>);
    }

    void collision_volume_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon, Vector8d& grad)
    {
        collision_constraint_grad_fd(Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul,
            impact_node, epsilon, collision_volume<double>, grad);
    }

    template double collision_volume<double>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, const ImpactNode, const double);

    template DScalar collision_volume<DScalar>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        DVector2 const&, DVector2 const&, DVector2 const&, DVector2 const&,
        const ImpactNode, const double);

    template void collision_volume_derivative<1>(const Eigen::MatrixX2d&,
        const Eigen::MatrixX2d&, const Eigen::MatrixX2i&, const EdgeEdgeImpact&,
        const int, const double, const int, Eigen::VectorXd&);

    template void collision_volume_derivative<Eigen::Dynamic>(
        const Eigen::MatrixX2d&, const Eigen::MatrixX2d&,
        const Eigen::MatrixX2i&, const EdgeEdgeImpact&, const int, const double,
        const int, Eigen::MatrixXd&);

} // namespace autodiff
} // namespace ccd
