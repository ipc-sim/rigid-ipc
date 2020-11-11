#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace physics {

    template <typename T> class Pose;
    template <typename T> using Poses = std::vector<Pose<T>>;

    /// @brief The position and rotation of an object.
    template <typename T> class Pose {
    public:
        Pose();
        Pose(
            const Eigen::VectorX3<T>& position,
            const Eigen::VectorX3<T>& rotation);
        Pose(const Eigen::VectorX6<T>& dof);
        Pose(const T& x, const T& y, const T& theta);
        // clang-format off
        Pose(const T& x, const T& y, const T& z,
             const T& theta_x, const T& theta_y, const T& theta_z);
        // clang-format on

        static Pose<T> Zero(int dim);

        static Poses<T> dofs_to_poses(const Eigen::VectorX<T>& dofs, int dim);
        static Eigen::VectorX<T> poses_to_dofs(const Poses<T>& poses);

        static int dim_to_ndof(const int dim) { return dim == 2 ? 3 : 6; }
        static int dim_to_pos_ndof(const int dim) { return dim; }
        static int dim_to_rot_ndof(const int dim)
        {
            return dim_to_ndof(dim) - dim_to_pos_ndof(dim);
        }
        int dim() const { return position.size(); }
        int pos_ndof() const { return position.size(); }
        int rot_ndof() const { return rotation.size(); }
        int ndof() const { return pos_ndof() + rot_ndof(); }

        Eigen::VectorX6<T> dof() const;

        /// @brief Replace a selected dof with the dof in other.
        /// @param R Rotation this rotations coordinates to the dof in
        ///          is_dof_selected.
        void select_dof(
            const Eigen::VectorX6b& is_dof_selected,
            const Pose<T>& other,
            const Eigen::MatrixXX3d& R);
        /// @brief Replace a selected dof with the dof in other.
        void select_dof(
            const Eigen::VectorX6b& is_dof_selected, const Pose<T>& other)
        {
            select_dof(
                is_dof_selected, other,
                Eigen::MatrixXX3d::Identity(rot_ndof(), rot_ndof()));
        }
        /// @brief Zero out the i-th dof if is_dof_zero(i) == true.
        void zero_dof(
            const Eigen::VectorX6b& is_dof_zero, const Eigen::MatrixXX3d& R);
        /// @brief Zero out the i-th dof if is_dof_zero(i) == true.
        void zero_dof(const Eigen::VectorX6b& is_dof_zero)
        {
            zero_dof(
                is_dof_zero,
                Eigen::MatrixXX3d::Identity(rot_ndof(), rot_ndof()));
        }

        Eigen::MatrixXX3<T> construct_rotation_matrix() const;

        static Pose<T>
        interpolate(const Pose<T>& pose0, const Pose<T>& pose1, T t);

        bool operator==(const Pose<T>& other) const;

        friend Pose<T> operator*(const Pose<T>& pose, const T& x)
        {
            return Pose<T>(pose.position * x, pose.rotation * x);
        }
        friend Pose<T> operator*(const T& x, const Pose<T>& pose)
        {
            return Pose<T>(x * pose.position, x * pose.rotation);
        }
        inline Pose<T>& operator*=(const T& x);
        Pose<T> operator/(const T& x) const;

        template <typename T1> Pose<T1> cast() const
        {
            return Pose<T1>(
                position.template cast<T1>(), rotation.template cast<T1>());
        }

        /// Position dof (either 2D or 3D)
        Eigen::VectorX3<T> position;
        /// Rotation dof (either 1D or 3D) expressed as a rotation vector
        Eigen::VectorX3<T> rotation;
    };

    template <typename T>
    Poses<T> interpolate(const Poses<T>& pose0, const Poses<T>& pose1, T t);
    template <typename T> Poses<T> operator*(const Poses<T>& poses, const T& x);
    /// @brief Cast poses element-wise.
    template <typename T1, typename T2>
    Poses<T2> cast(const Poses<T1>& poses_T1);

} // namespace physics

template <typename T>
Eigen::MatrixXX3<T> construct_rotation_matrix(const Eigen::VectorX3<T>& r);

template <typename T> Eigen::Matrix3<T> rotate_to_z(Eigen::Vector3<T> n);
template <typename T> Eigen::Matrix3<T> rotate_around_z(const T& theta);
template <typename T>
void decompose_to_z_screwing(
    const physics::Pose<T>& pose_t0,
    const physics::Pose<T>& pose_t1,
    Eigen::Matrix3<T>& R0,
    Eigen::Matrix3<T>& P,
    T& omega);

} // namespace ccd

#include "pose.tpp"
