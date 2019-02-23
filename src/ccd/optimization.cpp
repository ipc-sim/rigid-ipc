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
        // ϕ(x) = x^2 / s^3 - 3 * x^2 / s^2 + 3 * x / s
        double g = x_s * (3 + x_s * (-3 + x_s)); // Horner's method
        return 1 / g - 1;
    }

    double phi_log(double x) { return x > 0 ? -log(x) : INF_D; }

    double phi_hookean(double x, double s)
    {
        return x > 0 ? pow(log(x / s), 2) : INF_D;
    }

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
        std::function<double(Eigen::VectorXd)> f,
        std::function<Eigen::VectorXd(Eigen::VectorXd)> gradient,
        std::function<Eigen::MatrixXd(Eigen::VectorXd)> hessian,
        std::function<bool(Eigen::VectorXd)> constraint, double mu,
        double epsilon)
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
            assert(constraint(x));
        } while (g.squaredNorm() > epsilon); // Stop when the gradient is zero
        return x; // x is a local minimum now
    }

}

}
