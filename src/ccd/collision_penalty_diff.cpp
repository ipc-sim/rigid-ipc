#include <ccd/collision_penalty_diff.hpp>

#include <autodiff/finitediff.hpp>
// #include <autogen/collision_penalty.hpp>
#include <ccd/time_of_impact.hpp>

#include <algorithm>
#include <iostream>

namespace ccd {
namespace autodiff {

    // -----------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Penalty Derivatives
    // -----------------------------------------------------------------------------

    void compute_penalties_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::VectorXd& penalties)
    {
        penalties.resize(E.rows());
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            penalties(i) = edge_impact_map[i] < 0
                ? (10 * barrier_epsilon)
                : collision_penalty_refresh_toi(V, U, E,
                    ee_impacts[size_t(edge_impact_map[i])], int(i),
                    barrier_epsilon);
        }
    }

    void compute_penalties_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        Eigen::MatrixXd& penalties_grad)
    {
        penalties_grad.resize(V.size(), E.rows()); // n x m
        Eigen::VectorXd grad(V.size());
        EdgeEdgeImpact ee_impact;
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            grad.setZero();
            if (edge_impact_map[i] >= 0) {
                ee_impact = ee_impacts[size_t(edge_impact_map[i])];
                collision_penalty_grad(
                    V, U, E, ee_impact, int(i), barrier_epsilon, grad);
            }
            penalties_grad.col(i) = grad;
        }
    }

    void compute_penalties_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double barrier_epsilon,
        std::vector<Eigen::MatrixXd>& penalties_hessian)
    {
        penalties_hessian.clear();
        penalties_hessian.reserve(size_t(E.rows()));

        Eigen::MatrixXd hessian(V.size(), V.size());
        EdgeEdgeImpact ee_impact;
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            hessian.setZero();
            if (edge_impact_map[i] >= 0) {
                ee_impact = ee_impacts[size_t(edge_impact_map[i])];
                collision_penalty_hessian(
                    V, U, E, ee_impact, int(i), barrier_epsilon, hessian);
            }
            penalties_hessian.push_back(hessian);
        }
    }

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Penalties & Derivatives
    // -----------------------------------------------------------------------------
    double collision_penalty_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon)
    {
        ccd::autodiff::ImpactNode impact_node;
        // We need to figure out which one is the impact node
        Eigen::Vector2i e_ij = edges.row(edge_id);
        Eigen::Vector2i e_kl;

        if (edge_id == impact.impacted_edge_index) {
            // IJ is the impactED edge. The impacting node is K or L.
            e_kl = edges.row(impact.impacting_edge_index);
            impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vK
                                                       : ccd::autodiff::vL;

        } else if (edge_id == impact.impacting_edge_index) {
            // IJ is the impactING edge
            e_kl = edges.row(impact.impacted_edge_index);
            impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vI
                                                       : ccd::autodiff::vJ;
        } else {
            return 0.0;
        }

        // Then we get the position and velocity of the edge vertices
        Eigen::Vector2d Vi = vertices.row(e_ij(0));
        Eigen::Vector2d Vj = vertices.row(e_ij(1));
        Eigen::Vector2d Ui = displacements.row(e_ij(0));
        Eigen::Vector2d Uj = displacements.row(e_ij(1));

        Eigen::Vector2d Vk = vertices.row(e_kl(0));
        Eigen::Vector2d Vl = vertices.row(e_kl(1));
        Eigen::Vector2d Uk = displacements.row(e_kl(0));
        Eigen::Vector2d Ul = displacements.row(e_kl(1));

        return collision_penalty(
            Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, impact_node, barrier_epsilon);
    }

    void collision_penalty_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::VectorXd& gradient)
    {
        return collision_penalty_derivative(vertices, displacements, edges,
            impact, edge_id, barrier_epsilon, /*order=*/1, gradient);
    }

    void collision_penalty_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, Eigen::MatrixXd& hessian)
    {
        return collision_penalty_derivative(vertices, displacements, edges,
            impact, edge_id, barrier_epsilon, /*order=*/2, hessian);
    }

    template <int T>
    void collision_penalty_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id,
        const double barrier_epsilon, const int order,
        Eigen::Matrix<double, Eigen::Dynamic, T>& derivative)
    {
        ccd::autodiff::ImpactNode impact_node;
        assert(order == 1 || order == 2);
        assert(T == 1 || T == Eigen::Dynamic);
        assert(order != 2 || T == Eigen::Dynamic);

        // We need to figure out which one is the impact node
        Eigen::Vector2i e_ij = edges.row(edge_id);
        Eigen::Vector2i e_kl;

        if (edge_id == impact.impacted_edge_index) {
            // IJ is the impactED edge. The impacting node is K or L.
            e_kl = edges.row(impact.impacting_edge_index);
            impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vK
                                                       : ccd::autodiff::vL;

        } else if (edge_id == impact.impacting_edge_index) {
            // IJ is the impactING edge
            e_kl = edges.row(impact.impacted_edge_index);
            impact_node = impact.impacting_alpha < 0.5 ? ccd::autodiff::vI
                                                       : ccd::autodiff::vJ;
        } else {
            return;
        }

        // Then we get the position and velocity of the edge vertices
        Eigen::Vector2d Vi = vertices.row(e_ij(0));
        Eigen::Vector2d Vj = vertices.row(e_ij(1));
        Eigen::Vector2d Ui = displacements.row(e_ij(0));
        Eigen::Vector2d Uj = displacements.row(e_ij(1));

        Eigen::Vector2d Vk = vertices.row(e_kl(0));
        Eigen::Vector2d Vl = vertices.row(e_kl(1));
        Eigen::Vector2d Uk = displacements.row(e_kl(0));
        Eigen::Vector2d Ul = displacements.row(e_kl(1));

        // LOCAL gradient and hessian. Indices refer to the 4 vertices involded
        // in the collision
        DScalar v = collision_penalty_differentiable(
            Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, impact_node, barrier_epsilon);

        // Assemble into GLOBAL gradient and hessian
        auto num_vertices = vertices.rows();
        int nodes[4] = { e_ij(0), e_ij(1), e_kl(0), e_kl(1) };

        // Note: global gradient is sorted as x,x,x,...y,y,y
        // while local gradient is sorted as x,y,x,y,...,x,y
        if (order == 1) {
            derivative.resize(vertices.size(), 1);
            derivative.setZero();
            Vector8d el_grad = v.getGradient();

            for (int i = 0; i < 4; i++) {
                derivative(nodes[i], 0) = el_grad(2 * i);
                derivative(nodes[i] + num_vertices, 0) = el_grad(2 * i + 1);
            }

        } else {
            derivative.resize(vertices.size(), vertices.size());
            derivative.setZero();
            Matrix8d el_hessian = v.getHessian();
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++) {
                    derivative(nodes[i], nodes[j]) = el_hessian(2 * i, 2 * j);
                    derivative(nodes[i] + num_vertices, nodes[j])
                        = el_hessian(2 * i + 1, 2 * j);
                    derivative(nodes[i] + num_vertices, nodes[j] + num_vertices)
                        = el_hessian(2 * i + 1, 2 * j + 1);
                    derivative(nodes[i], nodes[j] + num_vertices)
                        = el_hessian(2 * i, 2 * j + 1);
                }
        }
    }

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Penalty Derivatives
    // -----------------------------------------------------------------------------
    DScalar collision_penalty_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DVector2 DUl = dvector(6, Ul);
        DScalar penalty(0.0);

        penalty = collision_penalty(
            Vi, Vj, Vk, Vl, DUi, DUj, DUk, DUl, impact_node, barrier_epsilon);
        return penalty;
    }

    template <typename T>
    T collision_penalty(const Eigen::Vector2d& Vi, const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& Vk, const Eigen::Vector2d& Vl,
        const Vector2T<T>& Ui, const Vector2T<T>& Uj, const Vector2T<T>& Uk,
        const Vector2T<T>& Ul, const ImpactNode impact_node,
        const double barrier_epsilon)
    {
        T toi, alpha;
        bool success;

        T t = T(1);                             // Final time for the time step
        T time_scale = t + 2 * barrier_epsilon; // 2x for safety

        const Vector2T<T> scaled_Ui = time_scale * Ui,
                          scaled_Uj = time_scale * Uj,
                          scaled_Uk = time_scale * Uk,
                          scaled_Ul = time_scale * Ul;

        switch (impact_node) {
        case vK:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vk, scaled_Ui, scaled_Uj, scaled_Uk, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vk, scaled_Ui, scaled_Uj, scaled_Uk, toi, alpha);
            break;
        case vL:
            success = compute_edge_vertex_time_of_impact(
                Vi, Vj, Vl, scaled_Ui, scaled_Uj, scaled_Ul, toi);
            success &= temporal_parameterization_to_spatial(
                Vi, Vj, Vl, scaled_Ui, scaled_Uj, scaled_Ul, toi, alpha);
            break;
        case vI:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vi, scaled_Uk, scaled_Ul, scaled_Ui, toi);
            alpha = T(0);
            break;
        case vJ:
            success = compute_edge_vertex_time_of_impact(
                Vk, Vl, Vj, scaled_Uk, scaled_Ul, scaled_Uj, toi);
            alpha = T(1);
            break;
        }

        if (!success) {
            return T(10 * barrier_epsilon);
        }

        return time_scale * toi - t;
    }

    void collision_penalty_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double barrier_epsilon, Vector8d& grad)
    {
        auto f = [&](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);
            Eigen::Vector2d ul = U.segment(6, 2);

            double penalty = collision_penalty(
                Vi, Vj, Vk, Vl, ui, uj, uk, ul, impact_node, barrier_epsilon);
            return penalty;
        };

        Eigen::VectorXd finite_diff;
        Vector8d x;
        x.segment(0, 2) = Ui;
        x.segment(2, 2) = Uj;
        x.segment(4, 2) = Uk;
        x.segment(6, 2) = Ul;
        ccd::finite_gradient(x, f, finite_diff);
        grad << finite_diff;
    }

    template double collision_penalty<double>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        Eigen::Vector2d const&, const ImpactNode, const double);

    template DScalar collision_penalty<DScalar>(Eigen::Vector2d const&,
        Eigen::Vector2d const&, Eigen::Vector2d const&, Eigen::Vector2d const&,
        DVector2 const&, DVector2 const&, DVector2 const&, DVector2 const&,
        const ImpactNode, const double);

    template void collision_penalty_derivative<1>(const Eigen::MatrixX2d&,
        const Eigen::MatrixX2d&, const Eigen::MatrixX2i&, const EdgeEdgeImpact&,
        const int, const double, const int, Eigen::VectorXd&);

    template void collision_penalty_derivative<Eigen::Dynamic>(
        const Eigen::MatrixX2d&, const Eigen::MatrixX2d&,
        const Eigen::MatrixX2i&, const EdgeEdgeImpact&, const int, const double,
        const int, Eigen::MatrixXd&);

} // namespace autodiff
} // namespace ccd
