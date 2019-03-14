#ifdef BUILD_WITH_OSQP

#include <opt/linearized_constraint_solver.hpp>
#include <opt/solver.hpp>

#include <Eigen/SparseCore>
#include <osqp.h>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <iostream>

#define INF_D (std::numeric_limits<double>::infinity())

namespace ccd {
namespace opt {

    bool solve_problem_with_linearized_constraints(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const double VOLUME_EPSILON, const DetectionMethod CCD_DETECTION_METHOD,
        const int max_iter, Eigen::MatrixX2d& Uopt)
    {
        assert(V.size() == U.size());
        assert(U.size() == Uopt.size());

        Eigen::MatrixXd U_flat = U;
        U_flat.resize(U_flat.size(), 1);

        // Load problem data
        const c_int NUM_VARIABLES = U_flat.size();

        // Quadratic Energy
        // || U - U0 ||^2 = (U - U0)^T * (U - U0) = U^T * U - 2 * U0^T * U + C
        // Quadratic matrix
        // P = 2 * I
        Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> P(
            NUM_VARIABLES, NUM_VARIABLES);
        P.setIdentity();
        P *= 2;
        // Quadratic linear term
        // q = - 2 * U0
        Eigen::VectorXd q = -2 * U_flat;

        Eigen::Matrix<c_float, 1, 1> x1 = (0.5 * U_flat.transpose() * P * U_flat
            + q.transpose() * U_flat + U_flat.transpose() * U_flat);
        c_float x2 = c_float((U_flat - U_flat).squaredNorm());
        assert(abs(x1[0] - x2) < 1e-8);

        // Linear constraints
        bool use_constraints = true;
        const c_int NUM_CONSTRAINTS
            = use_constraints ? E.rows() : NUM_VARIABLES;
        Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> A;
        Eigen::Matrix<c_float, Eigen::Dynamic, 1> l, u;
        if (use_constraints) {
            // Collision Volumes
            EdgeEdgeImpacts ee_impacts;
            Eigen::VectorXi edge_impact_map(E.rows());
            detect_collisions(
                V, U, E, CCD_DETECTION_METHOD, ee_impacts, edge_impact_map);
            Eigen::VectorXd volumes;
            Eigen::MatrixXd volume_gradient;
            compute_volumes_fixed_toi(
                V, U, E, ee_impacts, edge_impact_map, VOLUME_EPSILON, volumes);
            autodiff::compute_volumes_gradient(V, U, E, ee_impacts,
                edge_impact_map, VOLUME_EPSILON, volume_gradient);

            // V(U0) + ∇V(U0)(U - U0) ≥ 0 →
            // -V(U0) + ∇V(U0) * U0 ≤ ∇V(U0) * U ≤ ∞
            // Linear constraint matrix
            // A = ∇V(U0) ∈ R^(m × n)
            A = volume_gradient.transpose().sparseView();
            // Linear constraint lower bounds
            // ℓ = -V(U0) + ∇V(U0) * U0 ∈ R^m
            l = -volumes + volume_gradient.transpose() * U_flat;
            // Linear constraint upper bounds
            // u = ∞ ∈ R^m
            u = Eigen::Matrix<c_float, Eigen::Dynamic, 1>(NUM_CONSTRAINTS);
            u.setOnes();
            u *= INF_D;
        } else {
            A = Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int>(
                NUM_CONSTRAINTS, NUM_VARIABLES);
            A.setIdentity();
            // Linear constraint lower bounds
            l = Eigen::Matrix<c_float, Eigen::Dynamic, 1>(NUM_CONSTRAINTS);
            l.setOnes();
            l *= -INF_D;
            // Linear constraint upper bounds
            u = Eigen::Matrix<c_float, Eigen::Dynamic, 1>(NUM_CONSTRAINTS);
            u.setOnes();
            u *= INF_D;
        }

        // Populate data
        OSQPData data; // OSQPData
        data.n = NUM_VARIABLES;
        data.m = NUM_CONSTRAINTS;
        P.makeCompressed();
        data.P = csc_matrix(P.rows(), P.cols(), P.nonZeros(), P.valuePtr(),
            P.innerIndexPtr(), P.outerIndexPtr());
        data.q = q.data();
        A.makeCompressed();
        data.A = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(),
            A.innerIndexPtr(), A.outerIndexPtr());
        data.l = l.data();
        data.u = u.data();

        // Define Solver settings as default
        OSQPSettings settings;
        osqp_set_default_settings(&settings);
        settings.max_iter = max_iter;
        double epsilon = 1e-8;
        settings.eps_abs = epsilon;
        settings.eps_rel = epsilon;
        settings.eps_prim_inf = epsilon / 10;
        settings.eps_dual_inf = epsilon / 10;
        settings.polish = true;

        // Setup workspace
        OSQPWorkspace* work(osqp_setup(&data, &settings)); // Workspace

        // Solve Problem
        c_int exit_flag = osqp_solve(work);

        Eigen::Map<Eigen::VectorXd> x(work->solution->x, NUM_VARIABLES);

        Uopt.col(0)
            = Eigen::Map<Eigen::VectorXd>(work->solution->x, Uopt.rows());
        Uopt.col(1) = Eigen::Map<Eigen::VectorXd>(
            &(work->solution->x[Uopt.rows()]), Uopt.rows());

        // Check the solve was successful
        bool solve_successful = true;
        Eigen::VectorXd lhs = A * x;
        int i = 0;
        while (solve_successful && i < lhs.size()) {
            solve_successful &= l[i] <= (lhs[i] + 10 * epsilon)
                && (lhs[i] - 10 * epsilon) <= u[i];
            i++;
        }

        // Cleanup
        osqp_cleanup(work);
        return solve_successful;
    }

} // namespace opt
} // namespace ccd

#endif
