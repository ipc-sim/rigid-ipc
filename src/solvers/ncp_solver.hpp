#pragma once
#include <memory> // shared_ptr

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

#include <solvers/lcp_solver.hpp>
#include <solvers/optimization_solver.hpp>

#include <memory>

namespace ccd {
namespace opt {

    typedef std::function<void(const Eigen::VectorXd& x,
        const Eigen::VectorXd& alpha, const double gamma)>
        callback_intermediate_ncp;

    /**
     * @brief Method for updating candidate solution during NCP iterations
     *
     * Update: \f$ x_{i+1} = x_{i} + \Delta x \f$
     */
    enum class NcpUpdate {
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i\f$
        G_GRADIENT,
        /// \f$\Delta x = A^{-1} jac_x(x_i)^T \alpha_i + A^{-1}b - x_i\f$
        LINEARIZED
    };

    static const char* NcpUpdateNames[] = { "G_GRADIENT", "LINEARIZED" };

    class NCPSolver : public OptimizationSolver {
    public:
        NCPSolver();
        ~NCPSolver() override;

        OptimizationResults solve(OptimizationProblem& problem) override;
        bool solve_ncp(const Eigen::SparseMatrix<double>& hessian,
            const Eigen::VectorXd& b, OptimizationProblem& problem,
            Eigen::VectorXd& x_opt, Eigen::VectorXd& alpha_opt);
        void compute_linear_system(OptimizationProblem& problem);

        bool compute();
        void initialize();
        void solve_lcp(Eigen::VectorXd& delta_x);
        void move_to_unfeasible_domain(Eigen::VectorXd& delta_x, double& gamma);
        void update_candidate(Eigen::VectorXd& delta_x, double& gamma);

        Eigen::VectorXd Ainv(const Eigen::VectorXd& x);

        // ---------------------
        // Configuration
        // ---------------------
        bool keep_in_unfeasible;
        bool check_convergence;
        bool check_volume_increase;
        bool solve_for_active_cstr;
        double convergence_tolerance;
        NcpUpdate update_type;
        LCPSolver lcp_solver;

        // ---------------------
        // Optimization Specifics
        // ---------------------
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        OptimizationProblem* problem;
        //        callback_intermediate_ncp& callback;

        // -----------------------
        // Optimization Status
        // -----------------------
        std::shared_ptr<Eigen::SparseLU<Eigen::SparseMatrix<double>>> Asolver;
        Eigen::VectorXd g_xi;
        Eigen::VectorXi g_active;
        Eigen::SparseMatrix<double> jac_g_xi;

        // ----------------------
        // Optimization results
        // ----------------------
        Eigen::VectorXd xi;
        Eigen::VectorXd alpha_i;
    };

} // namespace opt
} // namespace ccd
