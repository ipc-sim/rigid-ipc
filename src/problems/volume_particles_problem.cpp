#include "volume_particles_problem.hpp"

#include <autogen/collision_volume.hpp>
#include <ccd/time_of_impact.hpp>

#include <logger.hpp>
#include <utils/flatten.hpp>

namespace ccd {

namespace opt {
    VolumeParticlesProblem::VolumeParticlesProblem(const std::string& name)
        : ParticlesProblem(name)
    {
    }

    void VolumeParticlesProblem::settings(const nlohmann::json& params)
    {
        constraint_.settings(params["volume_constraint"]);
        opt_solver_.settings(params["ncp_solver"]);
        opt_solver_.set_problem(*this);
        ParticlesProblem::settings(params["particles_problem"]);
    }

    Eigen::VectorXd VolumeParticlesProblem::eval_g(const Eigen::VectorXd& xk)
    {

        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        Eigen::VectorXd g_uk;
        constraint_.compute_constraints(uk, g_uk);
        return g_uk;
    }

    void VolumeParticlesProblem::eval_g(const Eigen::VectorXd& xk,
        Eigen::VectorXd& g_uk,
        Eigen::SparseMatrix<double>& g_uk_jacobian,
        Eigen::VectorXi& g_uk_active)
    {
        Eigen::MatrixXd uk = xk - vec_vertices_t0;
        ccd::unflatten(uk, 2);

        constraint_.compute_constraints(uk, g_uk, g_uk_jacobian, g_uk_active);
    }

} // namespace opt
} // namespace ccd
