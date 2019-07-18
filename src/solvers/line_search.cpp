// Functions for optimizing functions.
// Includes line search to find a step length to reduce a function.
#include "line_search.hpp"

#include <logger.hpp>

namespace ccd {
namespace opt {

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    bool line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        double& step_length,
        const double min_step_length)
    {
        return line_search(x, dir, f, Eigen::VectorXd::Zero(dir.size()),
            step_length, min_step_length, 0.0);
    }

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    bool line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        double& step_length,
        const double min_step_length,
        const double armijo_rule_coeff)
    {
        return constrained_line_search(
            x, dir, f, grad_fx, [](const Eigen::VectorXd&) { return true; },
            step_length, min_step_length, armijo_rule_coeff);
    }

    // Search along a search direction to find a scalar step_length in [0, 1]
    // such that f(x + step_length * dir) ≤ f(x).
    // TODO: Filter the dof that violate the constraints. These are the indices
    // i where ϕ([g(x)]_i) = ∞.
    bool constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const Eigen::VectorXd& grad_fx,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double& step_length,
        const double min_step_length,
        const double armijo_rule_coeff)
    {
        const double fx = f(x); // Function value we want to beat

        // Wolfe conditions:
        // Armijo rule
        const double wolfe1 = armijo_rule_coeff * dir.transpose() * grad_fx;
        auto armijo_rule = [&]() {
            if (armijo_rule_coeff != 0.0)
                return f(x + step_length * dir) <= fx + step_length * wolfe1;
            return f(x + step_length * dir) < fx;
        };
        // Curvature condition
        // double wolfe_c2 = 0.9;
        // const double wolfe2 = -wolfe_c2 * dir.transpose() * grad_fx;
        auto curvature_condition = [&]() {
            // return -dir.transpose * grad_f(x + step_length * dir) <= wolfe2;
            return true;
        };

        while (step_length >= min_step_length) {
            if (armijo_rule() && curvature_condition()
                && constraint(x + step_length * dir)) {
                return true;
            }
            step_length /= 2.0;
        }
        return false;
    }

    // Log samples along the search direction.
    void sample_search_direction(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f)
    {
        const auto sampling = Eigen::VectorXd::LinSpaced(51, -1, 1);
        for (int i = 0; i < sampling.size(); i++) {
            double step_length = sampling(i);
            spdlog::debug("method=line_search step_length={:g} obj={:.16g} "
                          "sqr_norm_grad={:.16g}",
                step_length, f(x + step_length * dir),
                grad_f(x + step_length * dir).squaredNorm());
        }
    }

} // namespace opt
} // namespace ccd
