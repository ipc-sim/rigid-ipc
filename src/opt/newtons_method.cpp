// Functions for optimizing functions.
// Includes Newton's method with and without constraints.
#include <cmath>
#include <iostream>

#include <Eigen/LU>
#include <Eigen/Sparse>

#include <opt/newtons_method.hpp>

namespace ccd {

namespace opt {

    // Search along a search direction to find a scalar gamma in [0, 1] such
    // that f(x + gamma * dir) ≤ f(x).
    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<bool(const Eigen::VectorXd&)>& constraint)
    {
        double gamma = 1.0;
        double fx = f(x);
        for (int i = 0; i <= 32; i++) {
            if (f(x + gamma * dir) < fx && constraint(x + gamma * dir))
                return gamma;
            gamma /= 2.0;
        }
        return 0.0;
    }

    // Perform a single step of Newton's Method to minimize a function f(x).
    double newtons_method_step(Eigen::VectorXd& x,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu, double epsilon)
    {
        Eigen::VectorXd g = gradient(x);
        Eigen::MatrixXd H = hessian(x);
        H.diagonal().array() += mu;
        Eigen::VectorXd delta_x = H.lu().solve(-g);
        double gamma = constrained_line_search(x, delta_x, f, constraint);
        if (gamma <= epsilon)
            return 0; // ToDo: Figure out a better way to return this case.
        x += gamma * delta_x; // Store the return value for x_{n+1} in x
        assert(constraint(x));
        return g.squaredNorm();
    }

    // Performa a Newton's Method to minimize a function f(x).
    void newtons_method(Eigen::VectorXd& x,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu, double epsilon, int max_iter)
    {
        int i = 0;
        while (i++ < max_iter
            && newtons_method_step(
                   x, f, gradient, hessian, constraint, mu, epsilon)
                > epsilon)
            ;
    }

}
}
