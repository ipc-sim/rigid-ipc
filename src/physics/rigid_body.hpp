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
            const Eigen::Vector3d& v,
            const Eigen::Vector3b& is_dof_fixed,
            const Eigen::Vector3d& x = Eigen::Vector3d::Zero(),
            const Eigen::Vector3d& x_prev = Eigen::Vector3d::Zero());

    public:
        // ------------------------------------------------------------------------
        // Factory Methods
        // ------------------------------------------------------------------------

        static RigidBody from_velocity(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::Vector3d& velocity,
            const Eigen::Vector3b& is_dof_fixed = Eigen::Vector3b::Zero())
        {
            return from_velocity(vertices, edges, Eigen::Vector3d::Zero(),
                velocity, is_dof_fixed);
        }

        static RigidBody from_velocity(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::Vector3d& position,
            const Eigen::Vector3d& velocity,
            const Eigen::Vector3b& is_dof_fixed = Eigen::Vector3b::Zero());

        static RigidBody from_displacement(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::Vector3d& displacement,
            const Eigen::Vector3b& is_dof_fixed = Eigen::Vector3b::Zero());

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

        Eigen::Matrix2d grad_theta(const double theta) const
        {
            Eigen::Matrix2d gradtheta_R;
            gradtheta_R << -sin(theta), -cos(theta), cos(theta), -sin(theta);
            return gradtheta_R;
        }
        // ------------------------------------------------------------------------
        // CCD Functions
        // ------------------------------------------------------------------------

        /// \brief: computes vertices position for given state
        /// returns the positions of all vertices in 'world space',
        /// taking into account the given body's position
        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> world_vertices(
            const Eigen::Matrix<T, 3, 1>& v) const;

        template <typename T>
        Eigen::Matrix<T, 2, 1> world_vertex(
            const Eigen::Matrix<T, 3, 1>& _position,
            const int vertex_idx) const;

        Eigen::MatrixXd world_vertices_gradient(
            const Eigen::Vector3d& velocity) const;
        std::vector<Eigen::Matrix3d> world_vertices_hessian(
            const Eigen::Vector3d& velocity) const;

        Eigen::MatrixXd world_vertices_gradient_finite(
            const Eigen::Vector3d& position) const;

        Eigen::MatrixXd world_vertices_gradient_exact(
            const Eigen::Vector3d& position) const;
        std::vector<Eigen::Matrix3d> world_vertices_hessian_exact(
            const Eigen::Vector3d& velocity) const;

        // ------------------------------------------------------------------------
        // Geometry
        // ------------------------------------------------------------------------
        Eigen::MatrixX2d vertices; ///< vertices positions in body space
        Eigen::MatrixX2i edges;    ///< vertices connectivity

        double mass;              ///< total mass (M) of the rigid body
        double moment_of_inertia; ///< moment of intertia (I) of the rigid body
        Eigen::Vector3b
            is_dof_fixed; ///< flag to indicate if dof is fixed (doesnt' change)
        Eigen::Matrix3d mass_matrix;
        Eigen::Matrix3d inv_mass_matrix;

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

    namespace RBDiff {
        static constexpr size_t NUM_VARS = 3;
        typedef Eigen::Vector3d Vector3d;
        typedef Eigen::Matrix3d Matrix3d;
        typedef DScalar1<double, Vector3d> DScalar1;
        typedef DScalar2<double, Vector3d, Matrix3d> DScalar2;
        typedef Eigen::Matrix<DScalar1, 3, 1> D1Vector3;
        typedef Eigen::Matrix<DScalar2, 3, 1> D2Vector3;
        typedef Eigen::Matrix<DScalar1, Eigen::Dynamic, Eigen::Dynamic>
            D1MatrixXd;
        typedef Eigen::Matrix<DScalar2, Eigen::Dynamic, Eigen::Dynamic>
            D2MatrixXd;
        typedef Eigen::Matrix<DScalar1, Eigen::Dynamic, 1> D1VectorXd;
        typedef Eigen::Matrix<DScalar2, Eigen::Dynamic, 1> D2VectorXd;

        inline void activate() { DiffScalarBase::setVariableCount(NUM_VARS); }

        inline D1Vector3 d1vars(const size_t i, const Vector3d& v)
        {
            return D1Vector3(DScalar1(i, v[0]), DScalar1(i + 1, v[1]),
                DScalar1(i + 2, v[2]));
        }

        inline D2Vector3 d2vars(const size_t i, const Vector3d& v)
        {
            return D2Vector3(DScalar2(i, v[0]), DScalar2(i + 1, v[1]),
                DScalar2(i + 2, v[2]));
        }

        inline Eigen::MatrixXd get_gradient(const D1VectorXd& x)
        {
            Eigen::MatrixXd grad(x.rows(), NUM_VARS);
            for (int i = 0; i < x.rows(); ++i) {
                grad.row(i) = x(i).getGradient();
            }
            return grad;
        }

        inline std::vector<Matrix3d> get_hessian(const D2VectorXd& x)
        {
            std::vector<Matrix3d> hess;
            hess.reserve(x.rows());

            for (int i = 0; i < x.rows(); ++i) {
                hess.push_back(x(i).getHessian());
            }
            return hess;
        }
    } // namespace RBDiff

} // namespace physics
} // namespace ccd

#include "rigid_body.tpp"

