/**
 *  We handle optimization problems of the form
 *      MIN     f(x)      x ∈ Rⁿ
 *
 *   s.t.       g_L ≤ g(x) ≤ g_U
 *              x_L ≤  x   ≤ x_U
 *
 */
#pragma once

#include <nlohmann/json.hpp>

#include <opt/optimization_problem.hpp>
#include <opt/optimization_results.hpp>

namespace ccd {
namespace opt {

    class OptimizationSolver {
    protected:
        Eigen::VectorXi free_dof; ///< @breif Indices of the free degrees.

    public:
        OptimizationSolver(const std::string& name);
        OptimizationSolver(const std::string& name, const int max_iterations);
        virtual ~OptimizationSolver();

        virtual void init_free_dof(Eigen::VectorXb is_dof_fixed);

        virtual OptimizationResults solve(OptimizationProblem& problem) = 0;

        /// \brief clear: reset all internal structures used
        virtual void clear()
        {
            throw NotImplementedError(
                "clear OptimizationSolver not implemented");
        }

        /// \brief clear: initializes the internal structures used to solve
        virtual void init(OptimizationProblem& /*problem*/)
        {
            throw NotImplementedError(
                "init OptimizationSolver not implemented");
        }

        /// \brief step_solve: takes one outer-step of the solver
        virtual OptimizationResults step_solve()
        {
            throw NotImplementedError(
                "step_solve OptimizationSolver not implemented");
        }
        virtual int num_outer_iterations()
        {
            throw NotImplementedError(
                "num_outer_iterations OptimizationSolver not implemented");
        }

        virtual void eval_f(
            const Eigen::MatrixXd& /*points*/, Eigen::VectorXd& /*fx*/)
        {
            throw NotImplementedError(
                "eval_f OptimizationSolver not implemented");
        }

        /// returns the kkt condition \grad f - grad_g^T \lambda = 0
        virtual Eigen::VectorXd get_grad_kkt() const
        {
            throw NotImplementedError(
                "get_grad_kkt OptimizationSolver not implemented");
        }

        virtual bool has_inner_solver() { return false; }

        virtual const OptimizationSolver& inner_solver()
        {
            throw NotImplementedError(
                "inner_solver OptimizationSolver not implemented");
        }

        int max_iterations; ///< @brief Maximum number of iteration.

        inline const std::string& name() const { return name_; }

        virtual void settings(const nlohmann::json& json)
        {
            max_iterations = json["max_iterations"].get<int>();
        }

        virtual nlohmann::json settings() const
        {
            nlohmann::json json;
            json["max_iterations"] = max_iterations;
            return json;
        }

    protected:
        std::string name_;
    };

} // namespace opt
} // namespace ccd
