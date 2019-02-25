/**
 * Functions for finding optimizing functions.
 * Includes Newton's method with and without constraints.
 */

#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include <Eigen/Core>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * g(x) = \frac{x^2}{s^3} - \frac{3x^2}{s^2} + \frac{3x}{s}
     * \end{align}
     * \begin{align}
     * \phi_{\text{spline}}(x) = \frac{1}{g(x)} - 1
     * \end{align}
     *
     * See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{spline}}(x)\f$.
     * @param s Denominator inside x which steepens the function
     *  \f$(s \in (0, 1])\f$
     */
    double phi_spline(double x, double s);
    double phi_spline_gradient(double x, double s);
    double phi_spline_hessian(double x, double s);

    /**
     * Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * \phi_{\log}(x) = -\log(x)
     * \end{align}
     *
     * See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\log}(x)\f$.
     */
    double phi_log(double x, double s);
    double phi_log_gradient(double x, double s);
    double phi_log_hessian(double x, double s);

    /**
     * Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * \phi_{\text{hookean}}(x) = \log^2\left(\frac{x}{s}\right)
     * \end{align}
     *
     * See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{hookean}}(x)\f$.
     * @param s Denominator inside x which steepens the function
     *  \f$(s \in (0, 1])\f$
     */
    double phi_hookean(double x, double s);
    double phi_hookean_gradient(double x, double s);
    double phi_hookean_hessian(double x, double s);

    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir, std::function<double(Eigen::VectorXd)> f,
        std::function<bool(Eigen::VectorXd)> constraint);

    Eigen::VectorXd newtons_method(const Eigen::VectorXd& x0,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu = 1e-5, double epsilon = 1e-12);

}
}

#endif
