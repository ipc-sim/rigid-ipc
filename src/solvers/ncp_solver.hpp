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
        const Eigen::VectorXd& alpha,
        const double gamma)>
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

    class NCPSolver : public virtual IStateOptimizationSolver {
    public:
        NCPSolver();
        NCPSolver(const std::string& name);
        ~NCPSolver() override;

        void set_problem(IVolumeProblem& problem) ;

        OptimizationResults solve() override;
        void init_solve() override;

//        OptimizationResults solve(IVolumeProblem& problem) ;
//        void init(IVolumeProblem& problem);

        void settings(const nlohmann::json& json) override;
        nlohmann::json settings() const override;
        const std::string& name() const override { return name_; }

        bool has_inner_solver() override { return false; }
        const IBarrierOptimizationSolver& inner_solver() override
        {
            throw NotImplementedError("inner_solver NCPSolver not implemented");
        }

        OptimizationResults step_solve() override;

        int num_outer_iterations() const override
        {
            return num_outer_iterations_;
        }
        Eigen::VectorXd get_grad_kkt() const override;

        // -----------------------------------------s

        bool solve_ncp(const Eigen::SparseMatrix<double>& hessian,
            const Eigen::VectorXd& b,
            IVolumeProblem& problem_ptr_,
            Eigen::VectorXd& x_opt,
            Eigen::VectorXd& alpha_opt);

        // ---------------------
        // Configuration
        // ---------------------
        bool do_line_search;
        bool solve_for_active_cstr;
        double convergence_tolerance;
        NcpUpdate update_type;
        LCPSolver lcp_solver;
        int max_iterations;

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
        Eigen::VectorXd lambda_i;

    protected:
        void compute_linear_system(IVolumeProblem& problem_ptr_);
        void compute_initial_solution();

        void solve_lcp(const Eigen::VectorXd& xi,
            const Eigen::VectorXd& gxi,
            const Eigen::SparseMatrix<double>& jac_gxi,
            const Eigen::VectorXi& gactive,
            Eigen::VectorXd& delta_x,
            Eigen::VectorXd& lambda_i) const;

        Eigen::VectorXd Ainv(const Eigen::VectorXd& x) const;

        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b;
        IVolumeProblem* problem_ptr_;

        int num_outer_iterations_;
        std::string name_;
    };
    void zero_out_fixed_dof(
        const Eigen::VectorXb& is_fixed, Eigen::SparseMatrix<double>& m);
} // namespace opt
} // namespace ccd
