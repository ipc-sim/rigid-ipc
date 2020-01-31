#pragma once

#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <utils/eigen_ext.hpp>

#include <physics/mass.hpp>
#include <physics/pose.hpp>

#include <autodiff/autodiff.h>

#include <iostream>

namespace ccd {
namespace physics {

    class RigidBody {

    protected:
        /**
         * @brief Create rigid body with center of mass at 0,0.
         *
         * @param vertices  Vertices of the rigid body in body space
         * @param faces     Vertices pairs defining the topology of the rigid
         *                  body
         */
        RigidBody(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& faces,
            const Eigen::MatrixXi& edges,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const double density,
            const Eigen::VectorXb& is_dof_fixed,
            const bool oriented);

    public:
        static RigidBody from_points(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& faces,
            const Eigen::MatrixXi& edges,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const double density,
            const Eigen::VectorXb& is_dof_fixed,
            const bool oriented);

        // Faceless version for convienence (useful for 2D)
        static RigidBody from_points(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& edges,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const double density,
            const Eigen::VectorXb& is_dof_fixed,
            const bool oriented)
        {
            return from_points(vertices, Eigen::MatrixXi(), edges, pose,
                velocity, density, is_dof_fixed, oriented);
        }

        enum Step { PREVIOUS_STEP = 0, CURRENT_STEP };

        // --------------------------------------------------------------------
        // State Functions
        // --------------------------------------------------------------------

        /// @brief: computes vertices position for current or previous state
        Eigen::MatrixXd world_vertices(const Step step = CURRENT_STEP) const
        {
            return world_vertices(step == PREVIOUS_STEP ? pose_prev : pose);
        }
        Eigen::MatrixXd world_vertices_t0() const
        {
            return world_vertices(PREVIOUS_STEP);
        }
        Eigen::MatrixXd world_vertices_t1() const
        {
            return world_vertices(CURRENT_STEP);
        }

        Eigen::MatrixXd world_velocities() const;

        // --------------------------------------------------------------------
        // CCD Functions
        // --------------------------------------------------------------------

        /// @brief: computes vertices position for given state
        /// returns the positions of all vertices in 'world space',
        /// taking into account the given body's position
        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> world_vertices(
            const Pose<T>& pose) const;
        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> world_vertices(
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& dof) const
        {
            return world_vertices(Pose<T>(dof));
        }

        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> world_vertex(
            const Pose<T>& pose, const int vertex_idx) const;
        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, 1> world_vertex(
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& dof,
            const int vertex_idx) const
        {
            return world_vertex<T>(Pose<T>(dof), vertex_idx);
        }

        template <>
        Eigen::VectorXd world_vertex(
            const Eigen::VectorXd& dof, const int vertex_idx) const
        {
            Pose<double> p(dof);
            Eigen::VectorXd v = world_vertex<double>(p, vertex_idx);
            // std::cout << "Vᵢ:" << vertices.row(vertex_idx) << std::endl
            //           << "p.position:" << p.position.transpose() << std::endl
            //           << "p.rotation:" << p.rotation.transpose() << std::endl
            //           << "wᵢ:" << v.transpose() << std::endl
            //           << "W(i):" <<
            //           world_vertices<double>(dof).row(vertex_idx)
            //           << std::endl
            //           << std::endl;
            return v;
        }

        Eigen::MatrixXd world_vertices_gradient(const Pose<double>& pose) const;

        Eigen::MatrixXd world_vertices_gradient_exact(
            const Pose<double>& pose) const;
        std::vector<Eigen::MatrixXd> world_vertices_hessian_exact(
            const Pose<double>& velocity) const;

        int dim() const { return vertices.cols(); }
        int ndof() const { return pose.ndof(); }
        int pos_ndof() const { return pose.pos_ndof(); }
        int rot_ndof() const { return pose.rot_ndof(); }

        // --------------------------------------------------------------------
        // Geometry
        // --------------------------------------------------------------------
        Eigen::MatrixXd vertices; ///< vertices positions in body space
        Eigen::MatrixXi faces;    ///< vertices connectivity
        Eigen::MatrixXi edges;    ///< vertices connectivity

        /// @brief total mass (M) of the rigid body
        double mass;
        /// @breif moment of intertia (I) of the rigid body
        Eigen::MatrixXd moment_of_inertia;
        /// @brief maximum distance from CM to a vertex
        double r_max;

        /// @brief Flag to indicate if dof is fixed (doesnt' change)
        Eigen::VectorXb is_dof_fixed;
        Eigen::MatrixXd mass_matrix;
        Eigen::MatrixXd inv_mass_matrix;

        bool is_oriented; ///< use edge orientation for normals

        // --------------------------------------------------------------------
        // State
        // --------------------------------------------------------------------
        /// @brief current timestep position and rotation of the center of mass
        Pose<double> pose;
        /// @brief previous timestep position and rotation of the center of mass
        Pose<double> pose_prev;

        /// @brief current timestep velocity of the center of mass
        Pose<double> velocity;
        /// @brief previous timestep velocity of the center of mass
        Pose<double> velocity_prev;
    };

} // namespace physics
} // namespace ccd

#include "rigid_body.tpp"
