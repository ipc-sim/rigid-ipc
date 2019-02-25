#include <Eigen/LU>
#include <ccd/optimization.hpp>
#include <cmath>
#include <iostream>

#define INF_D (std::numeric_limits<double>::infinity())

namespace ccd {

namespace opt {

    double phi_spline(double x, double s)
    {
        if (x <= 0)
            return INF_D;
        if (x >= s)
            return 0;
        double x_s = x / s;
        // g(x) = (x / s)^3 - 3 * (x / s)^2 + 3 * (x / s)
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        return 1 / g - 1;
    }

    double phi_spline_gradient(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        double x_s = x / s;
        // g(x) = (x / s)^3 - 3 * (x / s)^2 + 3 * (x / s)
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        // g'(x) = 3 * x^2 / s^3 - 6 * x / s^2 + 3 / s
        //       = (3 * (x / s)^2 - 6 * (x / s) + 3 ) / s
        double dg = (3 + x_s * (-6 + 3 * x_s)) / s; // Horner's method
        // ϕ'(x) = -g^-2 * g'
        return -1 / (g * g) * dg;
    }

    double phi_spline_hessian(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        double x_s = x / s;
        // g(x) = (x / s)^3 - 3 * (x / s)^2 + 3 * (x / s)
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        // g'(x) = 3 * x^2 / s^3 - 6 * x / s^2 + 3 / s
        //       = (3 * (x / s)^2 - 6 * (x / s) + 3 ) / s
        double dg = (3 + x_s * (-6 + 3 * x_s)) / s; // Horner's method
        // g''(x) = 6 * x / s^3 - 6 / s^2 = (6 * (x / s) - 6) / s^2
        double ddg = (6 * x_s - 6) / (s * s);
        // ϕ''(x) = 2g^-3 * (g')^2 + -g^-2 * g'' = (2 / g * (g')^2 - g'') / g^2
        return (2 / g * dg * dg - ddg) / (g * g);
    }

    double phi_log(double x, double s)
    {
        if (x <= 0)
            return INF_D;
        if (x >= s)
            return 0;
        return -log(x / s);
    }

    double phi_log_gradient(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        return -1 / x;
    }

    double phi_log_hessian(double x, double s)
    {
        if (x <= 0 || x >= s)
            return 0;
        return 1 / (x * x);
    }

    double phi_hookean(double x, double s)
    {
        // ϕ(x) = log(x / s) ^ 2
        return x > 0 ? pow(log(x / s), 2) : INF_D;
    }

    double phi_hookean_gradient(double x, double s)
    {
        if (x <= 0)
            return 0;
        // ϕ'(x) = 2 * log(x / s) / x
        return 2 * log(x / s) / x;
    };

    double phi_hookean_hessian(double x, double s)
    {
        if (x <= 0)
            return 0;
        // ϕ'(x) = 2 * log(x / s) / x = 2 * f(x) / g(x)
        double f = log(x / s);
        double df = 1 / x;
        double g = x;
        double dg = 1;
        return 2 * (g * df - f * dg) / (g * g); // Quotient rule
    };

    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir, std::function<double(Eigen::VectorXd)> f,
        std::function<bool(Eigen::VectorXd)> constraint)
    {
        double gamma = 1.0;
        double fx = f(x);
        for (int i = 0; i <= 32; i++) {
            if (f(x + gamma * dir) <= fx && constraint(x + gamma * dir))
                return gamma;
            gamma /= 2.0;
        }
        return 0.0;
    }

    Eigen::VectorXd newtons_method(const Eigen::VectorXd& x0,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu, double epsilon)
    {
        double gamma;
        Eigen::VectorXd x = x0, g, delta_x;
        Eigen::MatrixXd H,
            I = mu * Eigen::MatrixXd::Identity(x0.rows(), x0.rows());
        do {
            g = gradient(x);
            H = hessian(x) + I;
            delta_x = H.lu().solve(-g);
            gamma = constrained_line_search(x, delta_x, f, constraint);
            if (gamma <= epsilon)
                break;
            x += gamma * delta_x;
            bool test = constraint(x);
            assert(test);
        } while (g.squaredNorm() > epsilon); // Stop when the gradient is zero
        return x; // x is a local minimum now
    }

}

}
