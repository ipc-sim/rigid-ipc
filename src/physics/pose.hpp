#pragma once

#include <Eigen/Core>
#include <vector>

#include <utils/eigen_ext.hpp>

namespace ccd {
namespace physics {

    /// @brief The position and rotation of an object.
    template <typename T> class Pose {
    public:
        Pose();
        Pose(int dim);
        Pose(
            const Eigen::VectorX3<T>& position,
            const Eigen::VectorX3<T>& rotation);
        Pose(const Eigen::VectorX6<T>& dof);

        static std::vector<Pose<T>>
        dofs_to_poses(const Eigen::VectorX6<T>& dofs, int dim);
        static Eigen::VectorX6<T>
        poses_to_dofs(const std::vector<Pose<T>>& poses);

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

        Pose<T> operator+(Pose<T> other) const;
        inline Pose<T>& operator+=(Pose<T> other);
        Pose<T> operator-(Pose<T> other) const;
        inline Pose<T>& operator-=(Pose<T> other);
        Pose<T> operator/(T x) const;
        Pose<T> operator*(T x) const;
        static Pose<T> lerp_poses(Pose<T> pose0, Pose<T> pose1, T t);

        template <typename T1> Pose<T1> cast() const
        {
            return Pose<T1>(
                position.template cast<T1>(), rotation.template cast<T1>());
        }

        Eigen::VectorX3<T> position;
        // TODO: Use Quaternions
        // Eigen::Quaternion<T> rotation;
        Eigen::VectorX3<T> rotation; // Euler angles
    };

} // namespace physics
} // namespace ccd

#include "pose.tpp"
