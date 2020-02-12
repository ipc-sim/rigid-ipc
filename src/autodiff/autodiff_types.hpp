#ifndef AUTODIFF_TYPES_H
#define AUTODIFF_TYPES_H

#include <autodiff/autodiff.h>
#include <utils/eigen_ext.hpp>
#include <vector>

namespace ccd {

template <int N, int maxN = N> class AutodiffType {
public:
    typedef Eigen::Matrix<double, N, 1, Eigen::ColMajor, maxN, 1> VectorNd;
    typedef Eigen::Matrix<double, N, N, Eigen::ColMajor, maxN, maxN> MatrixNd;

    typedef DScalar1<double, VectorNd> DDouble1;
    typedef DScalar2<double, VectorNd, MatrixNd> DDouble2;

    typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, Eigen::Dynamic> D1MatrixXd;
    typedef Eigen::Matrix<DDouble2, Eigen::Dynamic, Eigen::Dynamic> D2MatrixXd;

    typedef Eigen::Matrix<DDouble1, Eigen::Dynamic, 1> D1VectorXd;
    typedef Eigen::Matrix<DDouble2, Eigen::Dynamic, 1> D2VectorXd;

    typedef Eigen::Matrix<DDouble1, 3, 1> D1Vector3d;
    typedef Eigen::Matrix<DDouble2, 3, 1> D2Vector3d;
    typedef Eigen::Matrix<DDouble1, 2, 1> D1Vector2d;
    typedef Eigen::Matrix<DDouble2, 2, 1> D2Vector2d;

    inline static void activate()
    {
        assert(N >= 0);
        DiffScalarBase::setVariableCount(N);
    }
    inline static void activate(size_t variableCount)
    {
        assert(N == Eigen::Dynamic);
        DiffScalarBase::setVariableCount(variableCount);
    }

    template <typename T>
    inline static Eigen::Matrix<T, Eigen::Dynamic, 1>
    dTvars(const size_t i, const Eigen::VectorXd& v)
    {
        Eigen::Matrix<T, Eigen::Dynamic, 1> vec;
        vec.resize(v.rows());
        for (int r = 0; r < v.rows(); r++) {
            vec[r] = T(i + r, v[r]);
        }
        return vec;
    }

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
        Eigen::MatrixXd grad(x.rows(), DiffScalarBase::getVariableCount());
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

} // namespace ccd
#endif
