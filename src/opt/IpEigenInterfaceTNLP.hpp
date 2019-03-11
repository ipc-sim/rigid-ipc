#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <IpTNLP.hpp>

namespace ccd {
/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    namespace ip {

        using namespace Ipopt;
        static const char* SolverReturnStrings[] = {
            "success",
            "maximum_iterations_exceeded",
            "maximum_cpu_time_exceeded",
            "stop_at_tiny_step",
            "stop_at_acceptable_point",
            "local_infeasibility",
            "user_requested_stop",
            "feasible_point_found",
            "diverging_iterates",
            "restoration_failed",
            "error_in_step_computation",
            "invalid_number_detected",
            "not_enough_degrees_of_freedom",
            "invalid_option",
            "insufficient_memory",
            "internal_error",
        };

        static const double NO_UPPER_BOUND = 2e19;
        static const double NO_LOWER_BOUND = -2e19;

        typedef std::function<double(const Eigen::VectorXd& x)> callback_f;
        typedef std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>
            callback_grad_f;
        typedef std::function<Eigen::VectorXd(const Eigen::VectorXd& x)>
            callback_g;
        typedef std::function<Eigen::MatrixXd(const Eigen::VectorXd& x)>
            callback_jac_g;
        typedef std::function<void(const Eigen::VectorXd& x,
            const double obj_value, const int iteration)>
            callback_intermediate;

        struct IpoptResult {
            Eigen::VectorXd x;
            double value;
            bool success;
            ApplicationReturnStatus app_status;
            SolverReturn status;
        };

        struct IpoptOptions {
            Eigen::VectorXd x;
            double value;
            bool success;
            ApplicationReturnStatus app_status;
            SolverReturn status;
        };

        IpoptResult minimize_ipopt(const callback_f f,
            const Eigen::VectorXd& x0, const callback_grad_f grad_f,
            const Eigen::VectorXd& x_lower,
            const Eigen::VectorXd& x_upper,
            const int num_constraints, const callback_g& g,
            const callback_jac_g& jac_g,
            const Eigen::VectorXd& g_lower,
            const Eigen::VectorXd& g_upper,
            const int verbosity = 5, const int max_iter = 3000,
            const callback_intermediate = nullptr);

        /**
         * @brief Class for interfacing IPOPT TNLP problem
         * with Eigen matrices.
         */
        class EigenInterfaceTNLP : public TNLP {

        public:
            EigenInterfaceTNLP(const callback_f f, const Eigen::VectorXd& x0,
                const callback_grad_f grad_f,
                const Eigen::VectorXd& x_lower = Eigen::VectorXd(),
                const Eigen::VectorXd& x_upper = Eigen::VectorXd(),
                const int num_constraints = 0, const callback_g& g = nullptr,
                const callback_jac_g& jac_g = nullptr,
                const Eigen::VectorXd& g_lower = Eigen::VectorXd(),
                const Eigen::VectorXd& g_upper = Eigen::VectorXd(),
                const callback_intermediate callback = nullptr);

            int num_vars;
            int num_constraints;
            Eigen::ArrayXd x0;
            Eigen::ArrayXd x_l, x_u, g_l, g_u;
            IpoptResult result;

            callback_f eval_f_;
            callback_grad_f eval_grad_f_;
            callback_g eval_g_;
            callback_jac_g eval_jac_g_;
            callback_intermediate intermediate_cb_;

            /**
             * n: (out) the number of variables in the problem (dimension of
             * $x$) m: (out) the number of constraints in the problem (dimension
             * of $ g(x)$) nnz_jac_g: (out), the number of nonzero entries in
             * the Jacobian of g. nnz_h_lag: (out), the number of nonzero
             * entries in the Hessian. index_style: (out), the numbering style
             * used for row/col entries in the sparse matrix format (C_STYLE:
             * 0-based, FORTRAN_STYLE: 1-based)
             */
            virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                Index& nnz_h_lag, IndexStyleEnum& index_style) override;

            /**
             * n: (in), the number of variables in the problem (dimension of $
             * x$). x_l: (out) the lower bounds $ x^L$ for $ x$. x_u: (out) the
             * upper bounds $ x^U$ for $ x$. m: (in), the number of constraints
             * in the problem (dimension of $ g(x)$). g_l: (out) the lower
             * bounds $ g^L$ for $ g(x)$. g_u: (out) the upper bounds $ g^U$ for
             * $ g(x)$. Note that EQUALITY constraints can be formulated in the
             * above formulation by setting the corresponding components of g_l
             * and g_u to the same value.
             * */
            virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                Index m, Number* g_l, Number* g_u) override;

            /**
             * n: (in), the number of variables in the problem (dimension of $
             * x$). init_x: (in), if true, this method must provide an initial
             * value for $ x$ x: (out), the initial values for the primal
             * variables, $ x$ init_z: (in), if true, this method must provide
             * an initial value for the bound multipliers $ z^L$ and $ z^U$ z_L:
             * (out), the initial values for the bound multipliers, $ z^L$. z_U:
             * (out), the initial values for the bound multipliers, $ z^U$. m:
             * (in), the number of constraints in the problem (dimension of $
             * g(x)$). init_lambda: (in), if true, this method must provide an
             * initial value for the constraint multipliers, $ \lambda$. lambda:
             * (out), the initial values for the constraint multipliers, $
             * \lambda$.
             */
            virtual bool get_starting_point(Index n, bool init_x, Number* x,
                bool init_z, Number* z_L, Number* z_U, Index m,
                bool init_lambda, Number* lambda) override;

            /**
             *
             * n: (in), the number of variables in the problem (dimension of $
             * x$). x: (in), the values for the primal variables, $ x$, at which
             * $ f(x)$ is to be evaluated. new_x: (in), false if any evaluation
             * method was previously called with the same values in x, true
             * otherwise. obj_value: (out) the value of the objective function
             * ($ f(x)$).
             * **/
            virtual bool eval_f(Index n, const Number* x, bool new_x,
                Number& obj_value) override;

            /**
             * Return the gradient of the objective function at the point $ x$.
             *
             * n: (in), the number of variables in the problem (dimension of $
             * x$). x: (in), the values for the primal variables, $ x$, at which
             * $ \nabla f(x)$ is to be evaluated. new_x: (in), false if any
             * evaluation method was previously called with the same values in
             * x, true otherwise. grad_f: (out) the array of values for the
             * gradient of the objective function ( $ \nabla f(x)$).
             */
            virtual bool eval_grad_f(
                Index n, const Number* x, bool new_x, Number* grad_f) override;

            /**
             * Return the value of the constraint function at the point $ x$.
             *
             * n: (in), the number of variables in the problem (dimension of $
             * x$). x: (in), the values for the primal variables, $ x$, at which
             * the constraint functions,  $ g(x)$, are to be evaluated. new_x:
             * (in), false if any evaluation method was previously called with
             * the same values in x, true otherwise. m: (in), the number of
             * constraints in the problem (dimension of $ g(x)$). g: (out) the
             * array of constraint function values, $ g(x)$.
             * */
            virtual bool eval_g(Index n, const Number* x, bool new_x, Index m,
                Number* g) override;

            /** Method to return:
             *   1) The structure of the jacobian (if "values" is NULL)
             *   2) The values of the jacobian (if "values" is not NULL)
             * n: (in), the number of variables in the problem (dimension of $
             * x$). x: (in), the values for the primal variables, $ x$, at which
             * the constraint Jacobian, $ \nabla g(x)^T$, is to be evaluated.
             * new_x: (in), false if any evaluation method was previously called
             * with the same values in x, true otherwise. m: (in), the number of
             * constraints in the problem (dimension of $ g(x)$). n_ele_jac:
             * (in), the number of nonzero elements in the Jacobian (dimension
             * of iRow, jCol, and values). iRow: (out), the row indices of
             * entries in the Jacobian of the constraints. jCol: (out), the
             * column indices of entries in the Jacobian of the constraints.
             * values: (out), the values of the entries in the Jacobian of the
             * constraints.
             */
            virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                Index m, Index nele_jac, Index* iRow, Index* jCol,
                Number* values) override;

            /** This method is called once per iteration, after the iteration
             *  summary output has been printed.  It provides the current
             *  information to the user to do with it anything she wants.  It
             *  also allows the user to ask for a premature termination of the
             *  optimization by returning FALSE, in which case Ipopt will
             *  terminate with a corresponding return status.  The basic
             *  information provided in the argument list has the quantities
             *  values printed in the iteration summary line.  If more
             *  information is required, a user can obtain it from the IpData
             *  and IpCalculatedQuantities objects.  However, note that the
             *  provided quantities are all for the problem that Ipopt sees,
             *  i.e., the quantities might be scaled, fixed variables might be
             *  sorted out, etc.  The status indicates things like whether the
             *  algorithm is in the restoration phase...  In the restoration
             *  phase, the dual variables are probably not not changing. */
            virtual bool intermediate_callback(AlgorithmMode mode, Index iter,
                Number obj_value, Number inf_pr, Number inf_du, Number mu,
                Number d_norm, Number regularization_size, Number alpha_du,
                Number alpha_pr, Index ls_trials, const IpoptData* ip_data,
                IpoptCalculatedQuantities* ip_cq) override;

            /** This method is called at the very end of the optimization.  It
             *  provides the final iterate to the user, so that it can be
             *  stored as the solution.  The status flag indicates the outcome
             *  of the optimization, where SolverReturn is defined in
             *  IpAlgTypes.hpp.  */
            virtual void finalize_solution(SolverReturn status, Index n,
                const Number* x, const Number* z_L, const Number* z_U, Index m,
                const Number* g, const Number* lambda, Number obj_value,
                const IpoptData* ip_data,
                IpoptCalculatedQuantities* ip_cq) override;
        };

    }
}
}
