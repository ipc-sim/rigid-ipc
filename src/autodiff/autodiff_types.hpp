#ifndef AUTODIFF_TYPES_H
#define AUTODIFF_TYPES_H

#include <autodiff/autodiff.h>
#include <vector>

namespace ccd {

template <int N> class AutodiffType {
public:
    typedef Eigen::Matrix<double, N, 1> VectorNd;
    typedef Eigen::Matrix<double, N, N> MatrixNd;

    typedef DScalar1<double, VectorNd> DDouble1;
    typedef DScalar2<double, VectorNd, MatrixNd> DDouble2;

    typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, Eigen::Dynamic> D1MatrixXd;
    typedef Eigen::Matrix<DDouble2, Eigen::Dynamic, Eigen::Dynamic> D2MatrixXd;

    typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, 1> D1VectorXd;
    typedef Eigen::Matrix<DDouble2, Eigen::Dynamic, 1> D2VectorXd;

    typedef Eigen::Matrix<DDouble1, 3, 1> D1Vector3d;
    typedef Eigen::Matrix<DDouble2, 3, 1> D2Vector3d;

    inline static void activate() { DiffScalarBase::setVariableCount(N); }

    inline static D1VectorXd d1vars(const size_t i, const Eigen::VectorXd& v)
    {
        D1VectorXd vec;
        vec.resize(v.rows());
        for (int r = 0; r < v.rows(); r++) {
            vec[r] = DDouble1(i + r, v[r]);
        }
        return vec;
    }
    inline static D2VectorXd d2vars(const size_t i, const Eigen::VectorXd& v)
    {
        D2VectorXd vec;
        vec.resize(v.rows());
        for (int r = 0; r < v.rows(); r++) {
            vec[r] = DDouble2(i + r, v[r]);
        }
        return vec;
    }
    inline static Eigen::MatrixXd get_gradient(const D1VectorXd& x)
    {
        Eigen::MatrixXd grad(x.rows(), N);
        for (int i = 0; i < x.rows(); ++i) {
            grad.row(i) = x(i).getGradient();
        }
        return grad;
    }

    inline static std::vector<MatrixNd> get_hessian(const D2VectorXd& x)
    {
        std::vector<MatrixNd> hess;
        hess.reserve(x.rows());

        for (int i = 0; i < x.rows(); ++i) {
            hess.push_back(x(i).getHessian());
        }
        return hess;
    }
};

///
/// \brief DistanceBarrierDiff the variables are 3 positions (x,y)
/// of an edge and a vertex.
///
typedef AutodiffType<6> DistanceBarrierDiff;

///
/// \brief TimeBarrierDiff the variables are 4 positions (x,y)
/// of two edges.
///
typedef AutodiffType<8> TimeBarrierDiff;

/// Used for time-barrier constraint.
/// TODO: should rewrite to use type defined above
static constexpr size_t N = 2 * 4;
typedef Eigen::Matrix<double, 8, 1> Vector8d;
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef DScalar2<double, Vector8d, Matrix8d> DScalar;
typedef DScalar::DVector2 DVector2;

static inline DVector2 dvector(
    size_t index, const Eigen::Matrix<double, 2, 1>& v)
{
    return DVector2(DScalar(index, v.x()), DScalar(index + 1, v.y()));
}

} // namespace ccd
#endif
