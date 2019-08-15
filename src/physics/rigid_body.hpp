#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <utils/eigen_ext.hpp>

#include <physics/center_of_mass.hpp>

#include <autodiff/autodiff.h>

namespace ccd {
namespace physics {

    class RigidBody {

    protected:
        ///
        /// \brief RigidBody:   Create rigid body with center of mass at 0,0
        ///
        /// \param vertices:    Vertices of the rigid body in body space
        /// \param edges:       Vertices pairs defining the topology
        ///                     of the rigid body
        /// \param v:           Velocity of the center of mass (v_x, v_y, omega)
        /// \param x:           Position and orientation of the center of mass
        /// \param x_prev:      Position and orientation of the center of mass
        /// of the previous step (x, y, theta)
        ///
        RigidBody(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXd& vertex_mass,
            const Eigen::Vector3b& is_dof_fixed,
            const bool oriented,
            const Eigen::Vector3d& velocity,
            const Eigen::Vector3d& position,
            const Eigen::Vector3d& position_prev);

    public:
        static RigidBody from_points(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXd& vertex_mass,
            const Eigen::Vector3b& is_dof_fixed,
            const bool oriented,
            const Eigen::Vector3d& position,
            const Eigen::Vector3d& velocity);

        enum Step { PREVIOUS_STEP = 0, CURRENT_STEP };

        // ------------------------------------------------------------------------
        // State Functions
        // ------------------------------------------------------------------------

        /// \brief: computes vertices position for current or previous state
        Eigen::MatrixXd world_vertices(const Step step = CURRENT_STEP) const;
        Eigen::MatrixXd world_vertices_t0() const
        {
            return world_vertices(PREVIOUS_STEP);
        }
        Eigen::MatrixXd world_vertices_t1() const
        {
            return world_vertices(CURRENT_STEP);
        }

        Eigen::MatrixXd world_velocities() const;

        Eigen::Matrix2d grad_theta(const double theta) const;

        // ------------------------------------------------------------------------
        // CCD Functions
        // ------------------------------------------------------------------------

        /// \brief: computes vertices position for given state
        /// returns the positions of all vertices in 'world space',
        /// taking into account the given body's position
        Eigen::MatrixXd world_vertices(const Eigen::Vector3d& v) const;

        template <typename T>
        Eigen::Matrix<T, 2, 1> world_vertex(
            const Eigen::Matrix<T, 3, 1>& position, const int vertex_idx) const;

        Eigen::MatrixXd world_vertices_gradient_exact(
            const Eigen::Vector3d& position) const;
        std::vector<Eigen::Matrix3d> world_vertices_hessian_exact(
            const Eigen::Vector3d& velocity) const;

        // ------------------------------------------------------------------------
        // Geometry
        // ------------------------------------------------------------------------
        Eigen::MatrixX2d vertices;       ///< vertices positions in body space
        Eigen::MatrixX2i edges;          ///< vertices connectivity
        Eigen::VectorXd per_vertex_mass; ///< vertices masses

        double mass;              ///< total mass (M) of the rigid body
        double moment_of_inertia; ///< moment of intertia (I) of the rigid body
        Eigen::Vector3b
            is_dof_fixed; ///< flag to indicate if dof is fixed (doesnt' change)
        Eigen::Matrix3d mass_matrix;
        Eigen::Matrix3d inv_mass_matrix;

        bool is_oriented; ///< use edge orientation for normals

        // ------------------------------------------------------------------------
        // State
        // ------------------------------------------------------------------------
        /// \brief current timestep velocity (v_x, v_y, v_theta)
        /// current timestep velocity of the center of mass
        Eigen::Vector3d velocity;

        /// \brief current timestep position (q_x, q_y, q_theta)
        /// position and orientation of the center of mass
        Eigen::Vector3d position;
        Eigen::Vector3d position_prev; ///> position of previous timestep
    };

} // namespace physics
} // namespace ccd

#include "rigid_body.tpp"
