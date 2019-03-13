#ifdef BUILD_WITH_OSQP

#include <opt/linearized_constraint_solve.hpp>

#include <Eigen/Sparse>
#include <osqp.h>

#include <ccd/collision_volume.hpp>
#include <ccd/collision_volume_diff.hpp>

#include <opt/displacements_opt.hpp>

#include <iostream>

#define INF_D (std::numeric_limits<double>::infinity())

namespace ccd {
namespace opt {

    void linearized_constraint_solve(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2d& U, const Eigen::MatrixX2i& E,
        const DetectionMethod CCD_DETECTION_METHOD, const double VOLUME_EPSILON)
    {
        assert(V.size() == U.size());

        Eigen::MatrixXd U_flat = U;
        U_flat.resize(U_flat.size(), 1);

        // Load problem data
        const c_int NUM_VARIABLES = U_flat.size();
        const c_int NUM_CONSTRAINTS = E.rows();

        // Quadratic Energy
        // || U - U0 ||^2 = (U - U0)^T * (U - U0) = U^T * U - 2 * U0^T * U + C
        // P = I; q = - 2 * U0
        // Quadratic matrix
        Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> P(
            NUM_VARIABLES, NUM_VARIABLES);
        P.setIdentity();
        // Quadratic linear term
        Eigen::VectorXd q = -2 * U_flat;

        // Collision Volumes
        EdgeEdgeImpacts ee_impacts;
        Eigen::VectorXi edge_impact_map(E.rows());
        detect_collisions(
            V, U, E, CCD_DETECTION_METHOD, ee_impacts, edge_impact_map);
        Eigen::VectorXd volumes;
        Eigen::MatrixXd volume_gradient;
        compute_volumes_fixed_toi(
            V, U, E, ee_impacts, edge_impact_map, VOLUME_EPSILON, volumes);
        autodiff::compute_volumes_gradient(V, U, E, ee_impacts, edge_impact_map,
            VOLUME_EPSILON, volume_gradient);

        // Linear constraints
        // V(U0) + ∇V(U0)(U - U0) ≥ 0 → -V(U0) + ∇V(U0) * U0 ≤ ∇V(U0) * U ≤ ∞
        // A = ∇V(U0) ∈ R^(m × n); ℓ = -V(U0) + ∇V(U0) * U0 ∈ R^m; u = ∞ ∈ R^m
        // Linear constraint matrix
        Eigen::SparseMatrix<c_float, Eigen::ColMajor, c_int> A(
            NUM_CONSTRAINTS, NUM_VARIABLES);
        A = volume_gradient.transpose().sparseView();
        A.makeCompressed();
        // Linear constraint lower bounds
        Eigen::VectorXd l = -volumes + volume_gradient.transpose() * U_flat;
        // Linear constraint upper bounds
        Eigen::Matrix<c_float, Eigen::Dynamic, 1> u(NUM_CONSTRAINTS);
        u.setOnes();
        u *= INF_D;

        // Populate data
        OSQPData data; // OSQPData
        data.n = NUM_VARIABLES;
        data.m = NUM_CONSTRAINTS;
        data.P = csc_matrix(P.rows(), P.cols(), P.nonZeros(), P.valuePtr(),
            P.innerIndexPtr(), P.outerIndexPtr());
        data.q = q.data();
        data.A = csc_matrix(A.rows(), A.cols(), A.nonZeros(), A.valuePtr(),
            A.innerIndexPtr(), A.outerIndexPtr());
        data.l = l.data();
        data.u = u.data();

        // Define Solver settings as default
        OSQPSettings settings;
        settings.alpha = 1.0; // Change alpha parameter
        osqp_set_default_settings(&settings);

        // Setup workspace
        OSQPWorkspace* work(osqp_setup(&data, &settings)); // Workspace

        // Solve Problem
        osqp_solve(work);

        Eigen::Map<Eigen::Matrix<c_float, Eigen::Dynamic, 1>> x(
            work->solution->x, NUM_VARIABLES);

        std::cout << "Solution:\n" << x << std::endl << std::endl;

        // Cleanup
        osqp_cleanup(work);
    }

} // namespace opt
} // namespace ccd

#endif
