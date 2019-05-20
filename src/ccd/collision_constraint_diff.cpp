// General code to compute the collision constraint and its derivative(s).
#include <ccd/collision_constraint_diff.hpp>

#include <Eigen/Sparse>

#include <autodiff/finitediff.hpp>
// #include <autogen/collision_volume.hpp>
#include <ccd/time_of_impact.hpp>

#include <iostream>

#include <profiler.hpp>
#ifdef PROFILE_FUNCTIONS
long number_of_constraint_calls = 0;
double time_spent_computing_constraint = 0;

long number_of_gradient_calls = 0;
double time_spent_computing_gradient = 0;

long number_of_hessian_calls = 0;
double time_spent_computing_hessian = 0;
#endif

namespace ccd {
namespace autodiff {

    ///
    /// \brief Determines which is the impacting node on an edge-edge impact
    /// \param[in] impact       : the edge-edge impact
    /// \param[in] edge_ij      : the edge we are evaluating
    /// \param[out] edge_kl     : the other edge involved in the impact
    /// \return                 : the impact node
    ///
    ImpactNode determine_impact_node(
        const EdgeEdgeImpact& impact, const int edge_ij, int& edge_kl)
    {
        if (edge_ij == impact.impacted_edge_index) {
            // IJ is the impactED edge. The impacting node is K or L.
            edge_kl = impact.impacting_edge_index;
            return impact.impacting_alpha < 0.5 ? ccd::autodiff::vK
                                                : ccd::autodiff::vL;

        } else if (edge_ij == impact.impacting_edge_index) {
            // IJ is the impactING edge
            edge_kl = impact.impacted_edge_index;
            return impact.impacting_alpha < 0.5 ? ccd::autodiff::vI
                                                : ccd::autodiff::vJ;
        } else {
            throw "Edge ID give is not part of the impact!";
        }
    }

    int get_constraint_index(
        const EdgeEdgeImpact& impact, const bool impacted, const int num_edges)
    {
        int e1 = impact.impacted_edge_index;
        int e2 = impact.impacting_edge_index;
        int p = impact.impacting_alpha > 0.5 ? 1 : 0;
        int q = impacted ? 0 : 1;

        // unravel index Q * P * E2 * e1 + Q * P * e2 + Q * p + q
        const int Q = 2, P = 2, E2 = num_edges;
        return Q * P * E2 * e1 + Q * P * e2 + Q * p + q;
    }

    int get_constraints_size(const int num_edges)
    {
        const int Q = 2, P = 2, E2 = num_edges, E1 = num_edges;
        return Q * P * E2 * E1;
    }

    // ------------------------------------------------------------------------
    // ALL IMPACTS GLOBAL Constraints Derivatives
    // ------------------------------------------------------------------------
    void compute_constraints_per_edge_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<double>& compute_constraint,
        const double default_value, Eigen::VectorXd& constraints)
    {
        constraints.resize(E.rows());
        constraints.setConstant(default_value);
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            if (edge_impact_map[i] >= 0) {
                constraints(i) = collision_constraint_refresh_toi(V, U, E,
                    ee_impacts[size_t(edge_impact_map[i])], int(i), epsilon,
                    compute_constraint);
            }
        }
    }

    void compute_constraints_gradient_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraints_grad)
    {
        constraints_grad.resize(V.size(), E.rows()); // n x m

        Eigen::SparseMatrix<double> grad;
        EdgeEdgeImpact ee_impact;
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            grad.setZero();
            if (edge_impact_map[i] >= 0) {
                ee_impact = ee_impacts[size_t(edge_impact_map[i])];
                collision_constraint_grad(V, U, E, ee_impact, int(i), epsilon,
                    compute_constraint, grad);
            }
            constraints_grad.col(i) = grad;
        }
    }

    void compute_constraints_hessian_per_edge(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& edge_impact_map, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        std::vector<Eigen::MatrixXd>& constraints_hessian)
    {
        constraints_hessian.clear();
        constraints_hessian.reserve(size_t(E.rows()));

        Eigen::SparseMatrix<double> hessian;
        EdgeEdgeImpact ee_impact;
        for (long i = 0; i < edge_impact_map.rows(); ++i) {
            hessian.setZero();
            if (edge_impact_map[i] >= 0) {
                ee_impact = ee_impacts[size_t(edge_impact_map[i])];
                collision_constraint_hessian(V, U, E, ee_impact, int(i),
                    epsilon, compute_constraint, hessian);
            }
            constraints_hessian.push_back(hessian);
        }
    }

    void compute_constraints_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints)
    {
#ifdef PROFILE_FUNCTIONS
        number_of_constraint_calls++;
        igl::Timer timer;
        timer.start();
#endif

        constraints.resize(int(2 * ee_impacts.size()));

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            constraints(2 * int(i))
                = collision_constraint_refresh_toi(V, U, E, ee_impact,
                    ee_impact.impacted_edge_index, epsilon, compute_constraint);

            constraints(2 * int(i) + 1) = collision_constraint_refresh_toi(V, U,
                E, ee_impact, ee_impact.impacting_edge_index, epsilon,
                compute_constraint);
        }

#ifdef PROFILE_FUNCTIONS
        timer.stop();
        time_spent_computing_constraint += timer.getElapsedTime();
#endif
    }

    void compute_constraints_dense_refresh_toi(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<double>& compute_constraint,
        Eigen::VectorXd& constraints)
    {
        const int num_edges = int(E.rows());
        const int num_constr = get_constraints_size(num_edges);
        constraints.resize(num_constr);
        constraints.setZero();

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            constraints[get_constraint_index(
                ee_impact, /*impacted=*/true, num_edges)]
                = collision_constraint_refresh_toi(V, U, E, ee_impact,
                    ee_impact.impacted_edge_index, epsilon, compute_constraint);

            constraints[get_constraint_index(
                ee_impact, /*impacted=*/false, num_edges)]
                = collision_constraint_refresh_toi(V, U, E, ee_impact,
                    ee_impact.impacting_edge_index, epsilon,
                    compute_constraint);
        }
    }

    void compute_constraints_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraint_grad)
    {
#ifdef PROFILE_FUNCTIONS
        number_of_gradient_calls++;
        igl::Timer timer;
        timer.start();
#endif

        constraint_grad.resize(V.size(), int(2 * ee_impacts.size()));

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];
            Eigen::SparseMatrix<double> grad;

            collision_constraint_grad(V, U, E, ee_impact,
                ee_impact.impacted_edge_index, epsilon, compute_constraint,
                grad);
            constraint_grad.col(2 * int(i)) = grad;

            collision_constraint_grad(V, U, E, ee_impact,
                ee_impact.impacting_edge_index, epsilon, compute_constraint,
                grad);
            constraint_grad.col(2 * int(i) + 1) = grad;
        }

#ifdef PROFILE_FUNCTIONS
        timer.stop();
        time_spent_computing_gradient += timer.getElapsedTime();
#endif
    }

    void compute_constraints_dense_gradient(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::MatrixXd& constraint_grad)
    {
        const int num_edges = int(E.rows());
        const int num_constr = get_constraints_size(num_edges);

        constraint_grad.resize(int(V.size()), num_constr);
        constraint_grad.setZero();

        Eigen::SparseMatrix<double> grad;
        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            auto& ee_impact = ee_impacts[i];

            collision_constraint_grad(V, U, E, ee_impact,
                ee_impact.impacted_edge_index, epsilon, compute_constraint,
                grad);
            constraint_grad.col(
                get_constraint_index(ee_impact, /*impacted=*/true, num_edges))
                = grad;

            collision_constraint_grad(V, U, E, ee_impact,
                ee_impact.impacting_edge_index, epsilon, compute_constraint,
                grad);
            constraint_grad.col(
                get_constraint_index(ee_impact, /*impacted=*/false, num_edges))
                = grad;
        }
    }

    void compute_constraints_hessian(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const EdgeEdgeImpacts& ee_impacts,
        const Eigen::VectorXi& /*edge_impact_map*/, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        std::vector<Eigen::SparseMatrix<double>>& constraint_hessian)
    {
#ifdef PROFILE_FUNCTIONS
        number_of_hessian_calls++;
        igl::Timer timer;
        timer.start();
#endif

        constraint_hessian.clear();
        constraint_hessian.reserve(2 * ee_impacts.size());

        Eigen::SparseMatrix<double> hessian;
        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            const EdgeEdgeImpact& ee_impact = ee_impacts[i];

            collision_constraint_hessian(V, U, E, ee_impact,
                ee_impact.impacted_edge_index, epsilon, compute_constraint,
                hessian);
            constraint_hessian.push_back(hessian);

            collision_constraint_hessian(V, U, E, ee_impact,
                ee_impact.impacting_edge_index, epsilon, compute_constraint,
                hessian);
            constraint_hessian.push_back(hessian);
        }

#ifdef PROFILE_FUNCTIONS
        timer.stop();
        time_spent_computing_hessian += timer.getElapsedTime();
#endif
    }

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT GLOBAL Constraints & Derivatives
    // -----------------------------------------------------------------------------
    double collision_constraint_refresh_toi(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_ij, const double epsilon,
        const constraint_func<double>& compute_constraint)
    {
        // We need to figure out which one is the impact node
        int edge_kl;
        ImpactNode impact_node;
        impact_node = determine_impact_node(impact, edge_ij, edge_kl);
        Eigen::Vector2i e_ij = edges.row(edge_ij);
        Eigen::Vector2i e_kl = edges.row(edge_kl);

        // Then we get the position and velocity of the edge vertices
        Eigen::Vector2d Vi = vertices.row(e_ij(0));
        Eigen::Vector2d Vj = vertices.row(e_ij(1));
        Eigen::Vector2d Ui = displacements.row(e_ij(0));
        Eigen::Vector2d Uj = displacements.row(e_ij(1));

        Eigen::Vector2d Vk = vertices.row(e_kl(0));
        Eigen::Vector2d Vl = vertices.row(e_kl(1));
        Eigen::Vector2d Uk = displacements.row(e_kl(0));
        Eigen::Vector2d Ul = displacements.row(e_kl(1));

        return compute_constraint(
            Vi, Vj, Vk, Vl, Ui, Uj, Uk, Ul, impact_node, epsilon);
    }

    void collision_constraint_grad(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& gradient)
    {
        return collision_constraint_derivative(vertices, displacements, edges,
            impact, edge_id, epsilon, /*order=*/1, compute_constraint,
            gradient);
    }

    void collision_constraint_hessian(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& hessian)
    {
        return collision_constraint_derivative(vertices, displacements, edges,
            impact, edge_id, epsilon, /*order=*/2, compute_constraint, hessian);
    }

    void collision_constraint_derivative(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2d& displacements, const Eigen::MatrixX2i& edges,
        const EdgeEdgeImpact& impact, const int edge_id, const double epsilon,
        const int order, const constraint_func<DScalar>& compute_constraint,
        Eigen::SparseMatrix<double>& derivative)
    {
        ccd::autodiff::ImpactNode impact_node;
        assert(order == 1 || order == 2);

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
        DScalar v = collision_constraint_differentiable(Vi, Vj, Vk, Vl, Ui, Uj,
            Uk, Ul, impact_node, epsilon, compute_constraint);

        // Assemble into GLOBAL gradient and hessian
        int num_vertices = int(vertices.rows());
        int nodes[4] = { e_ij(0), e_ij(1), e_kl(0), e_kl(1) };

        // Note: global gradient is sorted as x,x,x,...y,y,y
        // while local gradient is sorted as x,y,x,y,...,x,y
        if (order == 1) {
            derivative.resize(int(vertices.size()), 1);
            derivative.setZero();
            Vector8d el_grad = v.getGradient();

            for (int i = 0; i < 4; i++) {
                derivative.coeffRef(nodes[i], 0) = el_grad(2 * i);
                derivative.coeffRef(nodes[i] + num_vertices, 0)
                    = el_grad(2 * i + 1);
            }

        } else {
            derivative.resize(int(vertices.size()), int(vertices.size()));
            derivative.setZero();
            Matrix8d el_hessian = v.getHessian();
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++) {
                    derivative.coeffRef(nodes[i], nodes[j])
                        = el_hessian(2 * i, 2 * j);
                    derivative.coeffRef(nodes[i] + num_vertices, nodes[j])
                        = el_hessian(2 * i + 1, 2 * j);
                    derivative.coeffRef(
                        nodes[i] + num_vertices, nodes[j] + num_vertices)
                        = el_hessian(2 * i + 1, 2 * j + 1);
                    derivative.coeffRef(nodes[i], nodes[j] + num_vertices)
                        = el_hessian(2 * i, 2 * j + 1);
                }
        }
    }

    // -----------------------------------------------------------------------------
    // SINGLE IMPACT LOCAL Constraints Derivatives
    // -----------------------------------------------------------------------------
    DScalar collision_constraint_differentiable(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon,
        const constraint_func<DScalar>& compute_constraint)
    {
        // All definitions using DScalar must be done after setVariableCount
        DiffScalarBase::setVariableCount(8);

        DVector2 DUi = dvector(0, Ui);
        DVector2 DUj = dvector(2, Uj);
        DVector2 DUk = dvector(4, Uk);
        DVector2 DUl = dvector(6, Ul);

        return compute_constraint(
            Vi, Vj, Vk, Vl, DUi, DUj, DUk, DUl, impact_node, epsilon);
    }

    void collision_constraint_grad_fd(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj, const Eigen::Vector2d& Vk,
        const Eigen::Vector2d& Vl, const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj, const Eigen::Vector2d& Uk,
        const Eigen::Vector2d& Ul, const ImpactNode impact_node,
        const double epsilon, const constraint_func<double>& compute_constraint,
        Vector8d& grad)
    {
        auto f = [&](const Eigen::VectorXd& U) {
            Eigen::Vector2d ui = U.segment(0, 2);
            Eigen::Vector2d uj = U.segment(2, 2);
            Eigen::Vector2d uk = U.segment(4, 2);
            Eigen::Vector2d ul = U.segment(6, 2);

            return compute_constraint(
                Vi, Vj, Vk, Vl, ui, uj, uk, ul, impact_node, epsilon);
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

} // namespace autodiff
} // namespace ccd
