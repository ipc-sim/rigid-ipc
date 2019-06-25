#include <opt/particles_problem.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <nlohmann/json.hpp>

#include <physics/mass_matrix.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    ParticlesDisplProblem::ParticlesDisplProblem()
        : constraint(nullptr)
        , intermediate_callback(nullptr)
        , use_mass_matrix(true)
    {
    }

    ParticlesDisplProblem::~ParticlesDisplProblem() {}

    void ParticlesDisplProblem::initialize(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2i& E, const Eigen::MatrixX2d& U,
        CollisionConstraint& cstr)
    {
        vertices = V;
        edges = E;
        displacements = U;

        constraint = &cstr;
        constraint->initialize(V, E, U);
        initProblem();
    }

    void ParticlesDisplProblem::initProblem()
    {
        u_ = displacements;
        u_.resize(displacements.size(), 1);

        init_num_vars();
        x0.resize(num_vars);
        x_lower.resize(num_vars);
        x_upper.resize(num_vars);
        x_lower.setConstant(NO_LOWER_BOUND);
        x_upper.setConstant(NO_UPPER_BOUND);

        assert(is_dof_fixed.size() == num_vars);
        // is_dof_fixed.resize(num_vars);
        // is_dof_fixed.setZero();

        init_num_constraints();
        g_lower.resize(num_constraints);
        g_upper.resize(num_constraints);
        g_lower.setConstant(0.0);
        g_upper.setConstant(NO_UPPER_BOUND);

        is_collision_set_frozen = false;

        init_mass_matrix();
    }

    void ParticlesDisplProblem::init_num_vars()
    {
        num_vars = int(displacements.size());
    }

    void ParticlesDisplProblem::init_num_constraints()
    {
        // TODO: this is wrong, it will depend on constraint type
        // num_constraints = Eigen::Dynamic;
        num_constraints = constraint->number_of_constraints();
    }

    // Initalize the mass matrix based on the edge length of incident edges.
    void ParticlesDisplProblem::init_mass_matrix()
    {
        if (!use_mass_matrix) {
            mass_matrix
                = Eigen::SparseMatrix<double>(vertices.size(), vertices.size());
            mass_matrix.setIdentity();
            return;
        }

        Eigen::VectorXd vertex_masses;
        physics::mass_vector(vertices, edges, vertex_masses);
        // Repeat the mass vector to make a mass matrix per dof
        mass_matrix
            = Eigen::SparseDiagonal<double>(vertex_masses.replicate(2, 1));
    }

    ////////////////////////////////////////////////////////////////////////////
    // Objective function and its derivatives.

    double ParticlesDisplProblem::eval_f(const Eigen::VectorXd& x)
    {
        Eigen::VectorXd diff = x - u_;
        return 0.5 * diff.transpose() * mass_matrix * diff;
    }

    Eigen::VectorXd ParticlesDisplProblem::eval_grad_f(const Eigen::VectorXd& x)
    {
        return mass_matrix * (x - u_);
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_hessian_f(
        const Eigen::VectorXd& x)
    {
        return Eigen::MatrixXd(mass_matrix);
    }

    Eigen::SparseMatrix<double> ParticlesDisplProblem::eval_hessian_f_sparse(
        const Eigen::VectorXd& x)
    {
        return mass_matrix;
    }

    ////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////
    // Constraint function and its derivatives.

    Eigen::VectorXd ParticlesDisplProblem::eval_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        Eigen::VectorXd g_uk;
        if (!is_collision_set_frozen && constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        PROFILE(constraint->compute_constraints(Uk, g_uk),
            ProfiledPoint::COMPUTING_CONSTRAINTS);
        return g_uk;
    };

    void ParticlesDisplProblem::eval_g(const Eigen::VectorXd& x,
        Eigen::VectorXd& g_uk, Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        if (constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        constraint->compute_constraints(Uk, g_uk, g_uk_jacobian, g_uk_active);
    }

    Eigen::MatrixXd ParticlesDisplProblem::eval_jac_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        if (!is_collision_set_frozen && constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        Eigen::MatrixXd jac_gx;
        PROFILE(constraint->compute_constraints_jacobian(Uk, jac_gx),
            ProfiledPoint::COMPUTING_GRADIENT);

        return jac_gx;
    };

    void ParticlesDisplProblem::eval_jac_g(
        const Eigen::VectorXd& x, Eigen::SparseMatrix<double>& jac_gx)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        if (constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        constraint->compute_constraints_jacobian(Uk, jac_gx);
    };

    std::vector<Eigen::SparseMatrix<double>>
    ParticlesDisplProblem::eval_hessian_g(const Eigen::VectorXd& x)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        if (!is_collision_set_frozen && constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        PROFILE(constraint->compute_constraints_hessian(Uk, hess_gx);
                , ProfiledPoint::COMPUTING_HESSIAN);
        return hess_gx;
    }

    void ParticlesDisplProblem::eval_g_and_gdiff(const Eigen::VectorXd& x,
        Eigen::VectorXd& g_uk, Eigen::MatrixXd& g_uk_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian)
    {
        Eigen::MatrixXd Uk = x;
        Uk.resize(x.rows() / 2, 2);

        if (!is_collision_set_frozen && constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }
        std::vector<Eigen::SparseMatrix<double>> hess_gx;
        constraint->compute_constraints_and_derivatives(
            Uk, g_uk, g_uk_jacobian, g_uk_hessian);
    }

    ////////////////////////////////////////////////////////////////////////////

    bool ParticlesDisplProblem::eval_intermediate_callback(
        const Eigen::VectorXd& x)
    {
        if (intermediate_callback != nullptr) {
            Eigen::MatrixXd Uk = x;
            Uk.resize(x.rows() / 2, 2);
            return intermediate_callback(x, Uk);
        }

        return true;
    }

    void ParticlesDisplProblem::enable_line_search_mode(
        const Eigen::VectorXd& max_x)
    {
        Eigen::MatrixXd Uk = max_x;
        Uk.resize(max_x.rows() / 2, 2);

        if (!is_collision_set_frozen && constraint->update_collision_set) {
            constraint->detectCollisions(Uk);
        }

        this->is_collision_set_frozen = true;
    }

    void ParticlesDisplProblem::disable_line_search_mode()
    {
        this->is_collision_set_frozen = false;
    }

    // Check if a type, T, is differentiable (differentiable<T>())
    template <> bool differentiable<DScalar>() { return true; }
    template <> bool differentiable<double>() { return false; }

} // namespace opt
} // namespace ccd
