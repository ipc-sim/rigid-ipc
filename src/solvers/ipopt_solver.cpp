#ifdef BUILD_WITH_IPOPT
#include "ipopt_solver.hpp"

#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpTNLPAdapter.hpp"

namespace ccd {
namespace opt {

    IpoptSolver::IpoptSolver()
        : tolerance(1e-8)
        , print_level(0)
    {
        initialize();
    }

    IpoptSolver::IpoptSolver(
        double tolerance, int print_level, int max_iterations)
        : OptimizationSolver(max_iterations)
        , tolerance(tolerance)
        , print_level(print_level)
    {
        initialize();
    }

    void IpoptSolver::initialize()
    {
        app = IpoptApplicationFactory();
        app->Options()->SetStringValue("mu_strategy", "adaptive");
        app->Options()->SetStringValue(
            "hessian_approximation", "limited-memory");
        app->Options()->SetStringValue("derivative_test", "first-order");
    }

    OptimizationResults IpoptSolver::solve(OptimizationProblem& problem)
    {

        // update settings
        app->Options()->SetNumericValue("tol", tolerance);
        app->Options()->SetIntegerValue("print_level", print_level);
        app->Options()->SetIntegerValue("max_iter", max_iterations);

        ApplicationReturnStatus app_status;
        app_status = app->Initialize();

        SmartPtr<EigenInterfaceTNLP> nlp = new EigenInterfaceTNLP(problem);
        app_status = app->OptimizeTNLP(nlp);

        OptimizationResults result = nlp->result;
        result.success = (app_status == Solve_Succeeded);
        return result;
    }

    EigenInterfaceTNLP::EigenInterfaceTNLP(OptimizationProblem& problem)
        : problem(&problem)
    {

        result.x = problem.x0;
    }

    bool EigenInterfaceTNLP::get_nlp_info(Index& n,
        Index& m,
        Index& nnz_jac_g,
        Index& nnz_h_lag,
        IndexStyleEnum& index_style)
    {
        n = problem->num_vars;
        m = problem->num_constraints;

        // TODO: we need to know the constraints to fill this in
        // if we want to make it sparse
        nnz_jac_g = problem->num_constraints * problem->num_vars;

        // we are using hessian-approximation
        nnz_h_lag = 0;

        index_style = TNLP::C_STYLE;

        return true;
    }

    bool EigenInterfaceTNLP::get_bounds_info(Index n,
        Number* x_lower,
        Number* x_upper,
        Index m,
        Number* g_lower,
        Number* g_upper)
    {
        assert(n == problem->num_vars);
        assert(m == problem->num_constraints);

        Eigen::Map<Eigen::VectorXd>(&x_lower[0], problem->num_vars)
            = problem->x_lower;
        Eigen::Map<Eigen::VectorXd>(&x_upper[0], problem->num_vars)
            = problem->x_upper;

        Eigen::Map<Eigen::VectorXd>(&g_lower[0], problem->num_constraints)
            = problem->g_lower;
        Eigen::Map<Eigen::VectorXd>(&g_upper[0], problem->num_constraints)
            = problem->g_upper;

        return true;
    }

    bool EigenInterfaceTNLP::get_starting_point(Index n,
        bool init_x,
        double* x,
        bool init_z,
        double* /*z_L*/,
        double* /*z_U*/,
        Index /*m*/,
        bool init_lambda,
        double* /*lambda*/)
    {
        assert(n == problem->num_vars);

        // We assume we only have starting values for x
        assert(init_x == true);
        assert(init_z == false);
        assert(init_lambda == false);

        Eigen::Map<Eigen::VectorXd>(&x[0], n) = problem->x0;

        return true;
    }

    bool EigenInterfaceTNLP::eval_f(
        Index n, const double* x, bool /*new_x*/, double& obj_value)
    {
        assert(n == problem->num_vars);
        auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], n);
        obj_value = problem->eval_f(xk);

        return true;
    }

    bool EigenInterfaceTNLP::eval_grad_f(
        Index n, const Number* x, bool /*new_x*/, Number* grad_f_val)
    {
        assert(n == problem->num_vars);

        auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], n);
        Eigen::VectorXd grad = problem->eval_grad_f(xk);
        Eigen::Map<Eigen::MatrixXd>(grad_f_val, n, 1) = grad;

        return true;
    }

    bool EigenInterfaceTNLP::eval_g(
        Index n, const Number* x, bool /*new_x*/, Index m, Number* gval)
    {
        assert(n == problem->num_vars);
        assert(m == problem->num_constraints);

        auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], n);
        Eigen::VectorXd gk = problem->eval_g(xk);
        assert(gk.rows() == m);

        Eigen::Map<Eigen::MatrixXd>(gval, m, 1) = gk;

        return true;
    }

    bool EigenInterfaceTNLP::eval_jac_g(Index n,
        const Number* x,
        bool /*new_x*/,
        Index m,
        Index nele_jac,
        Index* iRow,
        Index* jCol,
        Number* values)
    {
        assert(n == problem->num_vars);
        assert(m == problem->num_constraints);
        assert(nele_jac == n * m); // dense

        if (values == nullptr) {
            // dense jacobian
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    iRow[i * n + j] = i;
                    jCol[i * n + j] = j;
                }
            }
        } else {
            auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], n);
            Eigen::MatrixXd jac_g_val;

            // the derivative of constraint g^{(i)} with respect to variable
            // x^{(j)} is placed in row i and column j.
            jac_g_val = problem->eval_jac_g(xk);

            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    values[i * n + j] = jac_g_val(i, j);
                }
            }
        }

        return true;
    }
    bool EigenInterfaceTNLP::intermediate_callback(AlgorithmMode mode,
        Index iter,
        Number obj_value,
        Number /*inf_pr*/,
        Number /*inf_du*/,
        Number /*mu*/,
        Number /*d_norm*/,
        Number /*regularization_size*/,
        Number /*alpha_du*/,
        Number alpha_pr,
        Index /*ls_trials*/,
        const IpoptData* ip_data,
        IpoptCalculatedQuantities* ip_cq)
    {
        if (mode == RegularMode) {
            Ipopt::TNLPAdapter* tnlp_adapter = nullptr;
            if (ip_cq != nullptr) {
                Ipopt::OrigIpoptNLP* orignlp;
                orignlp = dynamic_cast<OrigIpoptNLP*>(
                    GetRawPtr(ip_cq->GetIpoptNLP()));
                if (orignlp != nullptr)
                    tnlp_adapter
                        = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
            }
            Eigen::VectorXd primals(problem->num_vars);
            Eigen::VectorXd dualeqs(problem->num_constraints);
            tnlp_adapter->ResortX(*ip_data->curr()->x(), primals.data());
            tnlp_adapter->ResortG(*ip_data->curr()->y_c(),
                *ip_data->curr()->y_d(), dualeqs.data());

            problem->eval_intermediate_callback(primals);
        }
        return true;
    }

    void EigenInterfaceTNLP::finalize_solution(SolverReturn status,
        Index n,
        const Number* x,
        const Number* /*z_L*/,
        const Number* /*z_U*/,
        Index /*m*/,
        const Number* /*g*/,
        const Number* /*lambda*/,
        Number obj_value,
        const IpoptData* /*ip_data*/,
        IpoptCalculatedQuantities* /*ip_cq*/)
    {
        assert(n == problem->num_vars);
        result.x = Eigen::Map<const Eigen::VectorXd>(&x[0], n);
        result.minf = obj_value;
        std::cout << "solver_return: " << SolverReturnStrings[status]
                  << std::endl;
        std::cout << result.x.transpose() << std::endl;
    }

} // namespace opt
} // namespace ccd
#endif
