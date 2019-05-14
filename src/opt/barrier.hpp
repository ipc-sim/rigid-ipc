/**
 * Barrier functions that grow to infinity as x -> 0+. Includes gradient and
 * hessian functions, too. These barrier functions can be used to impose
 * inequlity constraints on a function.
 */

#pragma once

#include <Eigen/Core>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /// @brief Barriers for enforcing inequality constraints.
    enum BarrierType {
        SPLINE, ///< @brief A piecewise C2 function
        LOG,    ///< @brief Negative log piecewise C0 function
        HOOKEAN ///< @brief Log squared function
    };

    /**
     * @brief Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * g(x, s) = \frac{1}{s^3}x^3 - \frac{3}{s^2}x^2 + \frac{3}{s}x
     * \end{align}
     * \begin{align}
     * \phi_{\text{spline}}(x, s) = \frac{1}{g(x, s)} - 1
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{spline}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the barrier function at x.
     */
    template <typename T> T spline_barrier(T x, double s);

    /**
     * @brief Derivative of the spline_barrier function with respect to x.
     *
     * \begin{align}
     * g'(x, s) = \frac{3}{s^3}x^2 - \frac{6}{s^2}x + \frac{3}{s}
     * \end{align}
     * \begin{align}
     * \phi'_{\text{spline}}(x, s) = \frac{-1}{[g(x)]^2}g'(x)
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{spline}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the derivative of the barrier function at x.
     */
    double spline_barrier_gradient(double x, double s);

    /**
     * @brief Second derivative of the spline_barrier function with respect to
     * x.
     *
     * \begin{align}
     * g''(x, s) = \frac{6}{s^3}x - \frac{6}{s^2}
     * \end{align}
     * \begin{align}
     * \phi''_{\text{spline}}(x, s) = \frac{2}{[g(x)]^3}[g'(x)]^2 +
     * \frac{-1}{[g(x)]^2}g''(x) \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{spline}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the second derivative of the barrier function at x.
     */
    double spline_barrier_hessian(double x, double s);

    /**
     * @brief Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * \phi_{\log}(x, s) = -\ln\left(\frac{x}{s}\right)
     * \end{align}
     *
     *  where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\log}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the barrier function at x.
     */
    double log_barrier(double x, double s);

    /**
     * @brief Derivative of the log_barrier function with respect to x.
     *
     * \begin{align}
     * \phi'_{\log}(x, s) = \frac{-1}{x}
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\log}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the derivative of the barrier function at x.
     */
    double log_barrier_gradient(double x, double s);

    /**
     * @brief Second derivative of the log_barrier function with respect to x.
     *
     * \begin{align}
     * \phi''_{\log}(x, s) = \frac{1}{x^2}
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\log}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the second derivative of the barrier function at x.
     */
    double log_barrier_hessian(double x, double s);

    /**
     * @brief Function that grows to infinity as x approaches 0 from the right.
     *
     * \begin{align}
     * \phi_{\text{hookean}}(x) = \ln^2\left(\frac{x}{s}\right)
     * \end{align}
     *
     * See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     *
     * @param x The x value at which to evaluate \f$\phi_{\text{hookean}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the barrier function at x.
     */
    double hookean_barrier(double x, double s);

    /**
     * @brief Derivative of the hookean_barrier function with respect to x.
     *
     * \begin{align}
     * \phi'_{\text{hookean}}(x) = 2\ln\left(\frac{x}{s}\right) \frac{1}{x}
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\text{hookean}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the derivative of the barrier function at x.
     */
    double hookean_barrier_gradient(double x, double s);

    /**
     * @brief Second derivative of the hookean_barrier function with respect to
     * x.
     *
     * \begin{align}
     * \phi''_{\text{hookean}}(x) = \frac{2 - 2\ln\left(\frac{x}{s}\right)}{x^2}
     * \end{align}
     *
     * where \f$s\f$ is a constant. See <a
     * href="https://cims.nyu.edu/gcl/papers/LIM-2013-schueller.pdf">here</a>
     * for the original definition and a discussion.
     * @param x The x value at which to evaluate \f$\phi_{\text{hookean}}(x)\f$.
     * @param s Denominator inside x which steepens the function.
     * @return The value of the second derivative of the barrier function at x.
     */
    double hookean_barrier_hessian(double x, double s);

} // namespace opt
} // namespace ccd
