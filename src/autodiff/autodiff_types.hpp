#ifndef AUTODIFF_TYPES_H
#define AUTODIFF_TYPES_H

#include <autodiff/autodiff.h>

namespace ccd {

template <int N> class AutodiffType {
public:
    typedef Eigen::Matrix<double, N, 1> VectorNd;
    typedef Eigen::Matrix<double, N, N> MatrixNd;

    typedef DScalar1<double, VectorNd> DScalar1;
    typedef DScalar2<double, VectorNd, MatrixNd> DScalar2;

    typedef Eigen::Matrix<DScalar1, Eigen::Dynamic, Eigen::Dynamic> D1MatrixXd;
    typedef Eigen::Matrix<DScalar2, Eigen::Dynamic, Eigen::Dynamic> D2MatrixXd;

    typedef Eigen::Matrix<DScalar1, Eigen::Dynamic, 1> D1VectorXd;
    typedef Eigen::Matrix<DScalar2, Eigen::Dynamic, 1> D2VectorXd;

    inline static void activate() { DiffScalarBase::setVariableCount(N); }

    inline static D1VectorXd d1vars(const size_t i, const Eigen::VectorXd& v)
    {
        D1VectorXd vec;
        vec.resize(v.rows());
        for (int r = 0; r < v.rows(); r++) {
            vec[r] = DScalar1(i + r, v[r]);
        }
        return vec;
    }
    inline static D2VectorXd d2vars(const size_t i, const Eigen::VectorXd& v)
    {
        D2VectorXd vec;
        vec.resize(v.rows());
        for (int r = 0; r < v.rows(); r++) {
            vec[r] = DScalar2(i + r, v[r]);
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

typedef AutodiffType<6> DistanceBarrierDiff;

/// Used for time-barrier constraint.
/// TODO: should rewrite to use class defined above
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
