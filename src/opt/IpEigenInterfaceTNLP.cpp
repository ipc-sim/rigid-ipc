#ifdef BUILD_WITH_IPOPT

#include "IpEigenInterfaceTNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpTNLPAdapter.hpp"

namespace ccd {
namespace opt {
    namespace ip {
        IpoptResult minimize_ipopt(const callback_f f,
            const Eigen::VectorXd& x0, const callback_grad_f grad_f,
            const Eigen::VectorXd& x_lower, const Eigen::VectorXd& x_upper,
            const int num_constraints, const callback_g& g,
            const callback_jac_g& jac_g, const Eigen::VectorXd& g_lower,
            const Eigen::VectorXd& g_upper, const int verbosity,
            const int max_iter, const callback_intermediate callback)
        {

            SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
            app->Options()->SetNumericValue("tol", 1e-9);
            app->Options()->SetStringValue("mu_strategy", "adaptive");
            app->Options()->SetStringValue(
                "hessian_approximation", "limited-memory");
            app->Options()->SetStringValue("derivative_test", "first-order");
            app->Options()->SetIntegerValue("print_level", verbosity);
            app->Options()->SetIntegerValue("max_iter", max_iter);

            ApplicationReturnStatus app_status;
            app_status = app->Initialize();

            SmartPtr<EigenInterfaceTNLP> nlp
                = new EigenInterfaceTNLP(f, x0, grad_f, x_lower, x_upper,
                    num_constraints, g, jac_g, g_lower, g_upper, callback);
            app_status = app->OptimizeTNLP(nlp);

            IpoptResult result = nlp->result;
            result.app_status = app_status;
            result.success = (app_status == Solve_Succeeded);
            return result;
        }

        EigenInterfaceTNLP::EigenInterfaceTNLP(const callback_f f,
            const Eigen::VectorXd& x0, const callback_grad_f grad_f,
            const Eigen::VectorXd& x_lower, const Eigen::VectorXd& x_upper,
            const int num_constraints, const callback_g& g,
            const callback_jac_g& jac_g, const Eigen::VectorXd& g_lower,
            const Eigen::VectorXd& g_upper,
            const callback_intermediate callback)
            : num_vars(int(x0.rows()))
            , num_constraints(num_constraints)
            , x0(x0)
            , x_l(x_lower)
            , x_u(x_upper)
            , g_l(g_lower)
            , g_u(g_upper)
            , eval_f_(f)
            , eval_grad_f_(grad_f)
            , eval_g_(g)
            , eval_jac_g_(jac_g)
            , intermediate_cb_(callback)
        {
            assert(num_vars == x_l.rows() || x_l.size() == 0);
            assert(num_vars == x_u.rows() || x_u.size() == 0);
            assert(num_constraints == g_l.rows() || g_l.size() == 0);
            assert(num_constraints == g_u.rows() || g_u.size() == 0);
            assert((g == nullptr) == (jac_g == nullptr));
            assert(f != nullptr);
            assert(grad_f != nullptr);

            if (x_l.size() == 0) {
                x_l.resize(num_vars);
                x_l.setConstant(NO_LOWER_BOUND); // no-lower-bound
            }
            if (x_u.size() == 0) {
                x_u.resize(num_vars);
                x_u.setConstant(NO_UPPER_BOUND); // no-upper-bound
            }
            if (g_l.size() == 0) {
                g_l.resize(num_constraints);
                g_l.setConstant(NO_LOWER_BOUND); // no-lower-bound
            }
            if (g_u.size() == 0) {
                g_l.resize(num_constraints);
                g_l.setConstant(NO_UPPER_BOUND); // no-upper-bound
            }

            result.x = x0;
        }

        bool EigenInterfaceTNLP::get_nlp_info(Index& n, Index& m,
            Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
        {
            n = num_vars;
            m = num_constraints;

            // TODO: we need to know the constraints to fill this in
            // if we want to make it sparse
            nnz_jac_g = num_constraints * num_vars;

            // we are using hessian-approximation
            nnz_h_lag = 0;

            index_style = TNLP::C_STYLE;

            return true;
        }

        bool EigenInterfaceTNLP::get_bounds_info(Index n, Number* x_lower,
            Number* x_upper, Index m, Number* g_lower, Number* g_upper)
        {
            assert(n == num_vars);
            assert(m == num_constraints);

            Eigen::Map<Eigen::VectorXd>(&x_lower[0], num_vars) = this->x_l;
            Eigen::Map<Eigen::VectorXd>(&x_upper[0], num_vars) = this->x_u;

            Eigen::Map<Eigen::VectorXd>(&g_lower[0], num_constraints)
                = this->g_l;
            Eigen::Map<Eigen::VectorXd>(&g_upper[0], num_constraints)
                = this->g_u;

            return true;
        }

        bool EigenInterfaceTNLP::get_starting_point(Index n, bool init_x,
            double* x, bool init_z, double* /*z_L*/, double* /*z_U*/,
            Index /*m*/, bool init_lambda, double* /*lambda*/)
        {
            assert(n == num_vars);

            // We assume we only have starting values for x
            assert(init_x == true);
            assert(init_z == false);
            assert(init_lambda == false);

            Eigen::Map<Eigen::VectorXd>(&x[0], x0.rows()) = this->x0;

            return true;
        }

        bool EigenInterfaceTNLP::eval_f(
            Index n, const double* x, bool /*new_x*/, double& obj_value)
        {
            assert(n == num_vars);
            auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], num_vars);
            obj_value = eval_f_(xk);

            return true;
        }

        bool EigenInterfaceTNLP::eval_grad_f(
            Index n, const Number* x, bool /*new_x*/, Number* grad_f_val)
        {
            assert(n == num_vars);

            auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], num_vars);
            Eigen::VectorXd grad = eval_grad_f_(xk);
            Eigen::Map<Eigen::MatrixXd>(grad_f_val, n, 1) = grad;

            return true;
        }

        bool EigenInterfaceTNLP::eval_g(
            Index n, const Number* x, bool /*new_x*/, Index m, Number* gval)
        {
            assert(n == num_vars);
            assert(m == num_constraints);

            if (eval_g_) {
                auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], num_vars);
                Eigen::VectorXd gk = eval_g_(xk);
                assert(gk.rows() == m);

                Eigen::Map<Eigen::MatrixXd>(gval, m, 1) = gk;
            }

            return true;
        }

        bool EigenInterfaceTNLP::eval_jac_g(Index n, const Number* x,
            bool /*new_x*/, Index m, Index nele_jac, Index* iRow, Index* jCol,
            Number* values)
        {
            assert(n == num_vars);
            assert(m == num_constraints);
            assert(nele_jac == num_constraints * num_vars); // dense

            if (values == nullptr) {
                // dense jacobian
                for (int i = 0; i < num_constraints; ++i) {
                    for (int j = 0; j < num_vars; ++j) {
                        iRow[i * num_vars + j] = i;
                        jCol[i * num_vars + j] = j;
                    }
                }
            } else if (eval_jac_g_) {
                auto xk = Eigen::Map<const Eigen::VectorXd>(&x[0], num_vars);
                Eigen::MatrixXd jac_g_val;
                jac_g_val = eval_jac_g_(xk);

                for (int i = 0; i < num_constraints; ++i) {
                    for (int j = 0; j < num_vars; ++j) {
                        values[i * num_vars + j] = jac_g_val(i, j);
                    }
                }
            }

            return true;
        }

        bool EigenInterfaceTNLP::intermediate_callback(AlgorithmMode mode,
            Index iter, Number obj_value, Number /*inf_pr*/, Number /*inf_du*/,
            Number /*mu*/, Number /*d_norm*/, Number /*regularization_size*/,
            Number /*alpha_du*/, Number /*alpha_pr*/, Index /*ls_trials*/,
            const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
        {
            if (mode == RegularMode) {
                Ipopt::TNLPAdapter* tnlp_adapter = nullptr;
                if (ip_cq != nullptr) {
                    Ipopt::OrigIpoptNLP* orignlp;
                    orignlp = dynamic_cast<OrigIpoptNLP*>(
                        GetRawPtr(ip_cq->GetIpoptNLP()));
                    if (orignlp != nullptr)
                        tnlp_adapter = dynamic_cast<TNLPAdapter*>(
                            GetRawPtr(orignlp->nlp()));
                }
                Eigen::VectorXd primals(num_vars);
                Eigen::VectorXd dualeqs(num_constraints);
                tnlp_adapter->ResortX(*ip_data->curr()->x(), primals.data());
                tnlp_adapter->ResortG(*ip_data->curr()->y_c(),
                    *ip_data->curr()->y_d(), dualeqs.data());
                if (intermediate_cb_) {
                    intermediate_cb_(primals, obj_value, dualeqs, iter);
                }
            }
            return true;
        }

        void EigenInterfaceTNLP::finalize_solution(SolverReturn status,
            Index /*n*/, const Number* x, const Number* /*z_L*/,
            const Number* /*z_U*/, Index /*m*/, const Number* /*g*/,
            const Number* /*lambda*/, Number obj_value,
            const IpoptData* /*ip_data*/, IpoptCalculatedQuantities* /*ip_cq*/)
        {
            result.x = Eigen::Map<const Eigen::VectorXd>(&x[0], num_vars);
            result.value = obj_value;
            result.status = status;
            std::cout << "solver_return: " << SolverReturnStrings[status]
                      << std::endl;
            std::cout << result.x.transpose() << std::endl;
        }

    } // namespace ip

} // namespace opt

} // namespace ccd

#endif
