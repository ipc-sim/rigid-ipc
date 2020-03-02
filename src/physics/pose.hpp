#pragma once

#include <Eigen/Core>
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

        Eigen::MatrixXX3<T> construct_rotation_matrix() const;
        std::vector<Eigen::MatrixXX3<T>>
        construct_rotation_matrix_gradient() const;
        std::vector<std::vector<Eigen::MatrixXX3<T>>>
        construct_rotation_matrix_hessian() const;

        bool operator==(const Pose<T>& other) const;

        Pose<T> operator+(const Pose<T>& other) const;
        inline Pose<T>& operator+=(const Pose<T>& other);
        Pose<T> operator-(const Pose<T>& other) const;
        inline Pose<T>& operator-=(const Pose<T>& other);
        Pose<T> operator/(const T& x) const;
        Pose<T> operator*(const T& x) const;
        inline Pose<T>& operator*=(const T& x);
        static Pose<T>
        lerp(const Pose<T>& pose0, const Pose<T>& pose1, const T& t);

        template <typename T1> Pose<T1> cast() const
        {
            return Pose<T1>(
                position.template cast<T1>(), rotation.template cast<T1>());
        }

        /// Position dof (either 2D or 3D)
        Eigen::VectorX3<T> position;
        /// Rotation dof (either 1D or 3D) expressed in Euler angles
        Eigen::VectorX3<T> rotation;
        // TODO: Use Quaternions
        // Eigen::Quaternion<T> rotation;
    };

    template <typename T>
    Poses<T> operator+(const Poses<T>& poses0, const Poses<T>& poses1);
    template <typename T>
    Poses<T> operator-(const Poses<T>& poses0, const Poses<T>& poses1);
    template <typename T> Poses<T> operator*(const Poses<T>& poses, const T& x);
    /// @brief Cast poses element-wise.
    template <typename T1, typename T2>
    Poses<T2> cast(const Poses<T1>& poses_T1);

} // namespace physics
} // namespace ccd

#include "pose.tpp"
