#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <physics/center_of_mass.hpp>

#include <autodiff/autodiff.h>

namespace ccd {
namespace physics {

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

    class RigidBody {
    public:
        RigidBody(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::Vector2d& position,
            const Eigen::Vector3d& velocity);

        Eigen::MatrixXd world_displacements() const;

        template <typename T>
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> world_displacements(
            const Eigen::Matrix<T, 3, 1>& velocity) const;

        Eigen::MatrixX2d world_vertices() const;

        /// word_displacement_gradient: computes gradient for each
        /// displacement, xy.
        ///
        /// Returns: Matrix size (2n * 3)
        Eigen::MatrixXd world_displacements_gradient(
            const Eigen::Vector3d& velocity) const;

        /// word_displacement_gradient: computes hessian for each
        /// displacement, xy.
        ///
        /// Returns: Tensor size (2n * 3 * 3);
        ///          really a list of 3 * 3 matrices
        std::vector<Eigen::Matrix3d> world_displacements_hessian(
            const Eigen::Vector3d& velocity) const;

        Eigen::MatrixXd world_displacements_gradient_finite(
            const Eigen::Vector3d& velocity) const;

        std::vector<Eigen::MatrixX2d> world_displacements_gradient_exact(
            const Eigen::Vector3d& velocity) const;

        std::vector<std::vector<Eigen::MatrixX2d>>
        compute_world_displacements_hessian() const;

        /// \brief vertices as distances from the center of mass
        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;

        /// \brief position: position of the object in world space
        Eigen::Vector2d position;

        /// \brief velocity: linear and angular velicity
        Eigen::Vector3d velocity;

        ///
        /// \brief Centered: Factory method to create a RB with
        /// with position equal to current center of mass.
        /// \param vertices
        /// \param edges
        /// \param velocity
        /// \return
        ///
        static inline RigidBody Centered(const Eigen::MatrixXd& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::Vector3d& velocity)
        {
            Eigen::RowVector2d x = center_of_mass(vertices, edges);
            Eigen::MatrixX2d centered_vertices = vertices.rowwise() - x;
            return RigidBody(centered_vertices, edges, x, velocity);
        }
    };

} // namespace physics
} // namespace ccd
