#include "pose.hpp"

#include <Eigen/Geometry>
#include <tbb/tbb.h>

#include <utils/not_implemented_error.hpp>

namespace ccd {
namespace physics {

    template <typename T>
    Pose<T>::Pose()
        : position()
        , rotation()
    {
    }

    template <typename T>
    Pose<T>::Pose(
        const Eigen::VectorX3<T>& position, const Eigen::VectorX3<T>& rotation)
        : position(position)
        , rotation(rotation)
    {
    }

    template <typename T> Pose<T>::Pose(const Eigen::VectorX6<T>& dof)
    {
        if (dof.size() == dim_to_ndof(2)) {
            position = dof.head(dim_to_pos_ndof(2));
            rotation = dof.tail(dim_to_rot_ndof(2));
        } else if (dof.size() == dim_to_ndof(3)) {
            position = dof.head(dim_to_pos_ndof(3));
            rotation = dof.tail(dim_to_rot_ndof(3));
        } else {
            throw NotImplementedError("Unknown pose convertion for given ndof");
        }
    }

    template <typename T>
    Pose<T>::Pose(const T& x, const T& y, const T& theta)
        : Pose(Eigen::Vector2<T>(x, y), Eigen::Vector1<T>())
    {
        rotation << theta;
    }

    template <typename T>
    Pose<T>::Pose(
        const T& x,
        const T& y,
        const T& z,
        const T& theta_x,
        const T& theta_y,
        const T& theta_z)
        : Pose(
              Eigen::Vector3<T>(x, y, z),
              Eigen::Vector3<T>(theta_x, theta_y, theta_z))
    {
    }

    template <typename T> Pose<T> Pose<T>::Zero(int dim)
    {
        assert(dim == 2 || dim == 3);
        return Pose(
            Eigen::VectorX<T>::Zero(Pose<T>::dim_to_pos_ndof(dim)),
            Eigen::VectorX<T>::Zero(Pose<T>::dim_to_rot_ndof(dim)));
    }

    template <typename T>
    Poses<T> Pose<T>::dofs_to_poses(const Eigen::VectorX<T>& dofs, int dim)
    {
        int ndof = dim_to_ndof(dim);
        int num_poses = dofs.size() / ndof;
        assert(dofs.size() % ndof == 0);
        Poses<T> poses;
        poses.reserve(num_poses);
        for (int i = 0; i < num_poses; i++) {
            poses.emplace_back(dofs.segment(i * ndof, ndof));
        }
        return poses;
    }

    template <typename T>
    Eigen::VectorX<T> Pose<T>::poses_to_dofs(const Poses<T>& poses)
    {
        const int ndof = poses.size() ? poses[0].ndof() : 0;
        Eigen::VectorX<T> dofs(poses.size() * ndof);
        for (size_t i = 0; i < poses.size(); i++) {
            assert(poses[i].ndof() == ndof);
            dofs.segment(i * ndof, ndof) = poses[i].dof();
        }
        return dofs;
    }

    template <typename T> Eigen::VectorX6<T> Pose<T>::dof() const
    {
        Eigen::VectorX6<T> pose_dof(ndof());
        pose_dof.head(pos_ndof()) = position;
        pose_dof.tail(rot_ndof()) = rotation;
        return pose_dof;
    }

    template <typename T>
    Eigen::MatrixXX3<T> Pose<T>::construct_rotation_matrix() const
    {
        if (dim() == 2) {
            return Eigen::Rotation2D<T>(rotation(0)).toRotationMatrix();
        } else {
            typedef Eigen::Vector3<T> Vector3T;
            return (Eigen::AngleAxis<T>(rotation.z(), Vector3T::UnitZ())
                    * Eigen::AngleAxis<T>(rotation.y(), Vector3T::UnitY())
                    * Eigen::AngleAxis<T>(rotation.x(), Vector3T::UnitX()))
                .toRotationMatrix();
        }
    }

    template <typename T>
    std::vector<Eigen::MatrixXX3<T>>
    Pose<T>::construct_rotation_matrix_gradient() const
    {
        std::vector<Eigen::MatrixXX3<T>> grad_R(
            rot_ndof(), Eigen::MatrixXX3<T>(dim(), dim()));
        if (dim() == 2) {
            grad_R[0].row(0) << -sin(rotation(0)), -cos(rotation(0));
            grad_R[0].row(1) << cos(rotation(0)), -sin(rotation(0));
        } else {
            // Construct 3D rotation matricies
            Eigen::Matrix<T, 3, 3> Rx =
                Eigen::AngleAxis<T>(rotation.x(), Eigen::Vector3<T>::UnitX())
                    .toRotationMatrix();
            Eigen::Matrix<T, 3, 3> Ry =
                Eigen::AngleAxis<T>(rotation.y(), Eigen::Vector3<T>::UnitY())
                    .toRotationMatrix();
            Eigen::Matrix<T, 3, 3> Rz =
                Eigen::AngleAxis<T>(rotation.z(), Eigen::Vector3<T>::UnitZ())
                    .toRotationMatrix();

            // Construct gradient of each rotation matrix wrt its angle
            Eigen::Matrix<T, 3, 3> grad_Rx;
            // clang-format off
            grad_Rx <<
                0,          0,                  0,
                0, -sin(rotation.x()), -cos(rotation.x()),
                0,  cos(rotation.x()), -sin(rotation.x());
            Eigen::Matrix<T, 3, 3> grad_Ry;
            grad_Rx <<
                -sin(rotation.y()), 0,  cos(rotation.y()),
                         0,         0,          0,
                -cos(rotation.y()), 0, -sin(rotation.y());
            Eigen::Matrix<T, 3, 3> grad_Rz;
            grad_Rx <<
                -sin(rotation.z()), -cos(rotation.z()), 0,
                 cos(rotation.z()), -sin(rotation.z()), 0,
                         0,                 0,          0;
            // clang-format on
            grad_R[0] = Rz * Ry * grad_Rx;
            grad_R[1] = Rz * grad_Ry * Rx;
            grad_R[2] = grad_Rz * Ry * Rx;
        }
        return grad_R;
    }

    template <typename T>
    std::vector<std::vector<Eigen::MatrixXX3<T>>>
    Pose<T>::construct_rotation_matrix_hessian() const
    {
        std::vector<std::vector<Eigen::MatrixXX3<T>>> hess_R(
            rot_ndof(),
            std::vector<Eigen::MatrixXX3<T>>(
                rot_ndof(), Eigen::MatrixXX3<T>(dim(), dim())));

        if (dim() == 2) {
            hess_R[0][0] =
                -Eigen::Rotation2D<T>(rotation(0)).toRotationMatrix();
        } else {
            // Construct 3D rotation matricies
            Eigen::Matrix<T, 3, 3> Rx =
                Eigen::AngleAxis<T>(rotation.x(), Eigen::Vector3<T>::UnitX())
                    .toRotationMatrix();
            Eigen::Matrix<T, 3, 3> Ry =
                Eigen::AngleAxis<T>(rotation.y(), Eigen::Vector3<T>::UnitY())
                    .toRotationMatrix();
            Eigen::Matrix<T, 3, 3> Rz =
                Eigen::AngleAxis<T>(rotation.z(), Eigen::Vector3<T>::UnitZ())
                    .toRotationMatrix();

            // Construct gradient of each rotation matrix wrt its angle
            Eigen::Matrix<T, 3, 3> grad_Rx;
            // clang-format off
            grad_Rx <<
                0,          0,                  0,
                0, -sin(rotation.x()), -cos(rotation.x()),
                0,  cos(rotation.x()), -sin(rotation.x());
            Eigen::Matrix<T, 3, 3> grad_Ry;
            grad_Rx <<
                -sin(rotation.y()), 0,  cos(rotation.y()),
                         0,         0,          0,
                -cos(rotation.y()), 0, -sin(rotation.y());
            Eigen::Matrix<T, 3, 3> grad_Rz;
            grad_Rx <<
                -sin(rotation.z()), -cos(rotation.z()), 0,
                 cos(rotation.z()), -sin(rotation.z()), 0,
                         0,                 0,          0;
            // clang-format on
            hess_R[0][0] = Rz * Ry * -Rx;          // ∂R/∂x ∂R/∂x
            hess_R[0][1] = Rz * grad_Ry * grad_Rx; // ∂R/∂x ∂R/∂y
            hess_R[0][2] = grad_Rz * Ry * grad_Rx; // ∂R/∂x ∂R/∂z
            hess_R[1][0] = Rz * grad_Ry * grad_Rx; // ∂R/∂y ∂R/∂x
            hess_R[1][1] = Rz * -Ry * Rx;          // ∂R/∂y ∂R/∂y
            hess_R[1][2] = grad_Rz * grad_Ry * Rx; // ∂R/∂y ∂R/∂z
            hess_R[2][0] = grad_Rz * Ry * grad_Rx; // ∂R/∂z ∂R/∂x
            hess_R[2][1] = grad_Rz * grad_Ry * Rx; // ∂R/∂z ∂R/∂y
            hess_R[2][2] = -Rz * Ry * Rx;          // ∂R/∂z ∂R/∂z
        }
        return hess_R;
    }

    template <typename T> bool Pose<T>::operator==(const Pose<T>& other) const
    {
        return this->position == other.position
            && this->rotation == other.rotation;
    }

    template <typename T> Pose<T> Pose<T>::operator+(const Pose<T>& other) const
    {
        assert(dim() == other.dim());
        return Pose<T>(
            this->position + other.position, this->rotation + other.rotation);
    }

    template <typename T> Pose<T>& Pose<T>::operator+=(const Pose<T>& other)
    {
        assert(dim() == other.dim());
        this->position += other.position;
        this->rotation += other.rotation;
        return *this;
    }

    template <typename T> Pose<T> Pose<T>::operator-(const Pose<T>& other) const
    {
        assert(dim() == other.dim());
        return Pose<T>(
            this->position - other.position, this->rotation - other.rotation);
    }

    template <typename T> Pose<T>& Pose<T>::operator-=(const Pose<T>& other)
    {
        assert(dim() == other.dim());
        this->position -= other.position;
        this->rotation -= other.rotation;
        return *this;
    }

    template <typename T> Pose<T> Pose<T>::operator/(const T& x) const
    {
        return Pose<T>(this->position / x, this->rotation / x);
    }

    template <typename T> Pose<T> Pose<T>::operator*(const T& x) const
    {
        return Pose<T>(this->position * x, this->rotation * x);
    }

    template <typename T> Pose<T>& Pose<T>::operator*=(const T& x)
    {
        this->position *= x;
        this->rotation *= x;
        return *this;
    }

    template <typename T>
    Pose<T>
    Pose<T>::lerp_poses(const Pose<T>& pose0, const Pose<T>& pose1, const T& t)
    {
        assert(pose0.dim() == pose1.dim());
        return (pose1 - pose0) * t + pose0;
    }

    template <typename T>
    Poses<T> operator+(const Poses<T>& poses0, const Poses<T>& poses1)
    {
        Poses<T> sum = poses0;
        tbb::parallel_for(
            size_t(0), sum.size(), [&](size_t i) { sum[i] += poses1[i]; });
        return sum;
    }

    template <typename T>
    Poses<T> operator-(const Poses<T>& poses0, const Poses<T>& poses1)
    {
        Poses<T> diff = poses0;
        tbb::parallel_for(
            size_t(0), diff.size(), [&](size_t i) { diff[i] -= poses1[i]; });
        return diff;
    }

    template <typename T> Poses<T> operator*(const Poses<T>& poses, const T& x)
    {
        Poses<T> product = poses;
        tbb::parallel_for(
            size_t(0), product.size(), [&](size_t i) { product[i] *= x; });
        return product;
    }

    template <typename T1, typename T2>
    Poses<T2> cast(const Poses<T1>& poses_T1)
    {
        Poses<T2> poses_T2;
        poses_T2.reserve(poses_T1.size());
        for (int i = 0; i < poses_T1.size(); i++) {
            poses_T2.push_back(poses_T1[i].template cast<T2>());
        }
        return poses_T2;
    }

} // namespace physics
} // namespace ccd
