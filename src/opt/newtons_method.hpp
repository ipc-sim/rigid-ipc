/**
 * Functions for optimizing functions.
 * Includes Newton's method with and without constraints.
 */

#pragma once

#include <Eigen/Core>

namespace ccd {

/**
 * @namespace ccd::opt
 * @brief Functions for optimizing functions.
 */
namespace opt {

    /**
     * @brief Search along a search direction to find a scalar \f$\gamma \in [0,
     * 1]\f$ such that \f$f(x + \gamma \vec{dir}) \leq f(x)\f$.
     *
     * @param x Starting point for the line search.
     * @param dir Direction to search along.
     * @param f Function of x to minimize.
     * @param constraint Constrain on x such that constraint(x) must be true.
     * @return A scalar \f$\gamma\f$ that scales the search direction in an
     *      optimal way.
     */
    double constrained_line_search(const Eigen::VectorXd& x,
        const Eigen::VectorXd& dir, std::function<double(Eigen::VectorXd)> f,
        std::function<bool(Eigen::VectorXd)> constraint);

    /**
     * @brief Perform a single step of Newton's Method to minimize a function
     * \f$f(x)\f$.
     *
     * @param x The starting value for optimization and the place to store the
     * resulting updated value.
     * @param f A function \f$f: \mathbb{R}^n \rightarrow \mathbb{R}\f$ to
     * minimize.
     * @param gradient A function to compute the gradient of \f$f(x)\f$,
     * \f$\nabla f: \mathbb{R}^n \rightarrow \mathbb{R}^n\f$
     * @param hessian A function to compute the hessian of \f$f(x)\f$, \f$H(f):
     * \mathbb{R}^n \rightarrow \mathbb{R}^{n \times n}\f$
     * @param constraint A boolean function to check if the constraints are
     * satisfied.
     * @param mu A small value to add to the hessian diagonal to prevent it from
     * being singular.
     * @param epsilon A small value to check for close enough to zero.
     */
    double newtons_method_step(Eigen::VectorXd& x,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu = 1e-5, double epsilon = 1e-12);

    /**
     * Performa a Newton's Method to minimize a function \f$f(x)\f$.
     *
     * @param x The starting value for optimization and the place to store the
     * resulting minimum value.
     * @param f A function \f$f: \mathbb{R}^n \rightarrow \mathbb{R}\f$ to
     * minimize.
     * @param gradient A function to compute the gradient of \f$f(x)\f$,
     * \f$\nabla f: \mathbb{R}^n \rightarrow \mathbb{R}^n\f$
     * @param hessian A function to compute the hessian of \f$f(x)\f$, \f$H(f):
     * \mathbb{R}^n \rightarrow \mathbb{R}^{n \times n}\f$
     * @param constraint A boolean function to check if the constraints are
     * satisfied.
     * @param mu A small value to add to the hessian diagonal to prevent it from
     * being singular.
     * @param epsilon A small value to check for close enough to zero.
     */
    void newtons_method(Eigen::VectorXd& x,
        const std::function<double(const Eigen::VectorXd&)>& f,
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& gradient,
        const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& hessian,
        const std::function<bool(const Eigen::VectorXd&)>& constraint,
        double mu = 1e-5, double epsilon = 1e-12, int max_iter = 1000);

}
}
