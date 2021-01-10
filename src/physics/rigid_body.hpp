#pragma once

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include <physics/pose.hpp>
#include <utils/eigen_ext.hpp>

#include <BVH.hpp>
#include <utils/mesh_selector.hpp>

namespace ccd {
namespace physics {

    enum RigidBodyType { STATIC, KINEMATIC, DYNAMIC };

    NLOHMANN_JSON_SERIALIZE_ENUM(
        RigidBodyType,
        { { STATIC, "static" },
          { KINEMATIC, "kinematic" },
          { DYNAMIC, "dynamic" } });

    class RigidBody {
    protected:
        /**
         * @brief Create rigid body with center of mass at \f$\vec{0}\f$.
         *
         * @param vertices  Vertices of the rigid body in body space
         * @param faces     Vertices pairs defining the topology of the rigid
         *                  body
         */
        RigidBody(
            const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& edges,
            const Eigen::MatrixXi& faces,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const Pose<double>& force,
            const double density,
            const Eigen::VectorX6b& is_dof_fixed,
            const bool oriented,
            const int group_id,
            const RigidBodyType type = RigidBodyType::DYNAMIC);

    public:
        static RigidBody from_points(
            const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& edges,
            const Eigen::MatrixXi& faces,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const Pose<double>& force,
            const double density,
            const Eigen::VectorX6b& is_dof_fixed,
            const bool oriented,
            const int group_id,
            const RigidBodyType type = RigidBodyType::DYNAMIC);

        // Faceless version for convienence (useful for 2D)
        static RigidBody from_points(
            const Eigen::MatrixXd& vertices,
            const Eigen::MatrixXi& edges,
            const Pose<double>& pose,
            const Pose<double>& velocity,
            const Pose<double>& force,
            const double density,
            const Eigen::VectorX6b& is_dof_fixed,
            const bool oriented,
            const int group_id,
            const RigidBodyType type = RigidBodyType::DYNAMIC)
        {
            return from_points(
                vertices, edges, Eigen::MatrixXi(), pose, velocity, force,
                density, is_dof_fixed, oriented, group_id);
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

        /// @brief Computes vertices position for given state.
        /// @return The positions of all vertices in 'world space',
        ///         taking into account the given body's position.
        template <typename T>
        Eigen::MatrixX<T> world_vertices(
            const Eigen::MatrixXX3<T>& R, const Eigen::VectorX3<T>& p) const;
        template <typename T>
        Eigen::MatrixX<T> world_vertices(const Pose<T>& _pose) const
        {
            return world_vertices<T>(
                _pose.construct_rotation_matrix(), _pose.position);
        }
        template <typename T>
        Eigen::MatrixX<T> world_vertices(const Eigen::VectorX6<T>& dof) const
        {
            return world_vertices(Pose<T>(dof));
        }

        template <typename T>
        Eigen::VectorX3<T> world_vertex(
            const Eigen::MatrixXX3<T>& R,
            const Eigen::VectorX3<T>& p,
            const int vertex_idx) const;
        template <typename T>
        Eigen::VectorX3<T>
        world_vertex(const Pose<T>& _pose, const int vertex_idx) const
        {
            return world_vertex<T>(
                _pose.construct_rotation_matrix(), _pose.position, vertex_idx);
        }
        template <typename T>
        Eigen::VectorX3<T>
        world_vertex(const Eigen::VectorX6<T>& dof, const int vertex_idx) const
        {
            return world_vertex<T>(Pose<T>(dof), vertex_idx);
        }

        /// @warning Will not resize jac or hess, so make sure it is large
        /// enough.
        template <typename DScalar>
        Eigen::MatrixXd world_vertices_diff(
            const Pose<double>& pose,
            long rb_v0_i,
            Eigen::MatrixXd& V,
            Eigen::MatrixXd& jac,
            Eigen::MatrixXd& hess) const;

        double edge_length(int edge_id) const
        {
            return (vertices.row(edges(edge_id, 1))
                    - vertices.row(edges(edge_id, 0)))
                .norm();
        }

        long num_vertices() const { return vertices.rows(); }
        long num_edges() const { return edges.rows(); }
        long num_faces() const { return faces.rows(); }
        int dim() const { return vertices.cols(); }
        int ndof() const { return pose.ndof(); }
        int pos_ndof() const { return pose.pos_ndof(); }
        int rot_ndof() const { return pose.rot_ndof(); }

        void compute_bounding_box(
            const physics::Pose<double>& pose_t0,
            const physics::Pose<double>& pose_t1,
            Eigen::VectorX3d& box_min,
            Eigen::VectorX3d& box_max) const;
        void compute_bounding_box(
            const physics::Pose<double>& pose,
            Eigen::VectorX3d& box_min,
            Eigen::VectorX3d& box_max) const
        {
            return compute_bounding_box(pose, pose, box_min, box_max);
        }

        // --------------------------------------------------------------------
        // Properties
        // --------------------------------------------------------------------

        /// @brief Group id of this body
        int group_id;

        /// @brief Dyanmic type of rigid body
        RigidBodyType type;

        // --------------------------------------------------------------------
        // Geometry
        // --------------------------------------------------------------------
        Eigen::MatrixXd vertices; ///< Vertices positions in body space
        Eigen::MatrixXi edges;    ///< Vertices connectivity
        Eigen::MatrixXi faces;    ///< Vertices connectivity

        double average_edge_length; ///< Average edge length

        /// @brief total mass (M) of the rigid body
        double mass;
        /// @brief moment of inertia measured with respect to the principal axes
        Eigen::VectorX3d moment_of_inertia;
        /// @brief rotation from the principal axes to the input orientation
        Eigen::MatrixXX3d R0;
        /// @brief maximum distance from CM to a vertex
        double r_max;
        /// @brief the mass matrix of the rigid body
        Eigen::DiagonalMatrixX6d mass_matrix;

        /// @brief Flag to indicate if dof is fixed (doesnt' change)
        Eigen::VectorX6b is_dof_fixed;

        /// @brief Use edge orientation for normal in 2D restitution
        bool is_oriented;

        /// @brief Local space BVH initalized at construction
        BVH::BVH bvh;
        MeshSelector mesh_selector;

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

        Eigen::MatrixXX3d Qdot;

        /// @brief external force acting on the body
        Pose<double> force;
    };

} // namespace physics
} // namespace ccd

#include "rigid_body.tpp"
