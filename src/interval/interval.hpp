// An interval object.
#pragma once

#include <utils/eigen_ext.hpp>

#include <filib/fi_lib.h>

#include <string>

namespace ipc::rigid {

class Interval {
public:
    Interval();
    Interval(double x);
    Interval(double x, double y);

    static Interval empty();

    /* ------------------------------------------------------------------- */
    /* --- IO (input/output)                                           --- */
    /* ------------------------------------------------------------------- */

    friend std::istream& operator>>(std::istream& is, Interval& a);
    friend std::ostream& operator<<(std::ostream& os, const Interval& a);

    /* ------------------------------------------------------------------- */
    /* --- interval arithmetic (basic operations)                      --- */
    /* ------------------------------------------------------------------- */

    friend Interval operator+(const Interval& a, const Interval& b);
    friend Interval operator+(const Interval& a, const double b);
    friend Interval operator+(const double a, const Interval& b);
    friend Interval operator+(const Interval& a);

    Interval& operator+=(const Interval& rhs);
    Interval& operator+=(const double& rhs);

    friend Interval operator-(const Interval& a, const Interval& b);
    friend Interval operator-(const Interval& a, const double b);
    friend Interval operator-(const double a, const Interval& b);
    friend Interval operator-(const Interval& a);

    Interval& operator-=(const Interval& rhs);
    Interval& operator-=(const double rhs);

    friend Interval operator*(const Interval& a, const Interval& b);
    friend Interval operator*(const Interval& a, double b);
    friend Interval operator*(double a, const Interval& b);

    Interval& operator*=(const Interval& rhs);
    Interval& operator*=(const double rhs);

    friend Interval operator/(const Interval& a, const Interval& b);
    friend Interval operator/(const Interval& a, const double b);
    friend Interval operator/(const double a, const Interval& b);

    Interval& operator/=(const Interval& rhs);
    Interval& operator/=(const double& rhs);

    /* ------------------------------------------------------------------- */
    /* --- Interval arithmetic (logical operations)                    --- */
    /* ------------------------------------------------------------------- */

    friend bool operator==(const Interval& a, const Interval& b);
    friend bool operator==(const Interval& a, const double b);

    friend bool operator!=(const Interval& a, const Interval& b);

    friend bool operator<(const Interval& a, const Interval& b);
    friend bool operator<(const double a, const Interval& b);
    friend bool operator<=(const Interval& a, const Interval& b);
    friend bool operator<=(const double a, const Interval& b);

    friend bool operator>(const Interval& a, const Interval& b);
    friend bool operator>(const Interval& a, const double b);
    friend bool operator>=(const Interval& a, const Interval& b);
    friend bool operator>=(const Interval& a, const double b);

    friend Interval operator|(const Interval& a, const Interval& b);
    friend Interval operator&(const Interval& a, const Interval& b);

    Interval& operator|=(const Interval& b);
    Interval& operator&=(const Interval& b);

    friend bool in(double a, const Interval& b);
    friend bool in(const Interval& a, const Interval& b);

    bool is_empty() const;
    bool zero_in() const;

    bool overlaps(const Interval& b) const;

    /* ------------------------------------------------------------------- */
    /* --- utilities, mid, diam, ...                                   --- */
    /* ------------------------------------------------------------------- */

    double mid() const;
    std::pair<Interval, Interval> bisect() const;
    friend bool disjoint(const Interval& a, const Interval& b);
    double diam() const;
    double drel() const;
    double width() const { return diam(); }
    Interval blow(double eps) const;

    // min max function, same as what BOOST does
    friend Interval max(const Interval& x, const Interval& y);
    friend Interval max(const Interval& x, double y);
    friend Interval max(double x, const Interval& y);

    friend Interval min(const Interval& x, const Interval& y);
    friend Interval min(const Interval& x, double y);
    friend Interval min(double x, const Interval& y);

    /* ------------------------------------------------------------------- */
    /* --- Interval arithmetic (elementary functions)                  --- */
    /* ------------------------------------------------------------------- */

    friend Interval exp(const Interval& a);
    friend Interval expm(const Interval& a);
    friend Interval sinh(const Interval& a);
    friend Interval cosh(const Interval& a);
    friend Interval coth(const Interval& a);
    friend Interval tanh(const Interval& a);
    friend Interval log(const Interval& a);
    friend Interval ln(const Interval& a);
    friend Interval lg1p(const Interval& a);
    friend Interval sqrt(const Interval& a);
    friend Interval sqr(const Interval& a);
    friend Interval asnh(const Interval& a);
    friend Interval asinh(const Interval& a);
    friend Interval acsh(const Interval& a);
    friend Interval acosh(const Interval& a);
    friend Interval acth(const Interval& a);
    friend Interval acoth(const Interval& a);
    friend Interval atnh(const Interval& a);
    friend Interval atanh(const Interval& a);
    friend Interval asin(const Interval& a);
    friend Interval acos(const Interval& a);
    friend Interval acot(const Interval& a);
    friend Interval atan(const Interval& a);
    friend Interval sin(const Interval& a);
    friend Interval cos(const Interval& a);
    friend Interval cot(const Interval& a);
    friend Interval tan(const Interval& a);
    friend Interval exp2(const Interval& a);
    friend Interval ex10(const Interval& a);
    friend Interval log2(const Interval& a);
    friend Interval lg10(const Interval& a);
    friend Interval erf(const Interval& a);
    friend Interval erfc(const Interval& a);

    friend Interval abs(const Interval& a);

    const double& lower() const { return INF; }
    const double& upper() const { return SUP; }

    double& lower() { return INF; }
    double& upper() { return SUP; }

protected:
    // Helper functions to convert from pure FILib interval to my Interval
    Interval(interval i)
    {
        INF = i.INF;
        SUP = i.SUP;
    }

    // Helper functions to convert from my Interval to pure FILib interval
    operator interval() const { return { INF, SUP }; }

    double INF, SUP;
};

template <typename Derived>
inline Eigen::VectorXd width(const Eigen::MatrixBase<Derived>& x)
{
    Eigen::VectorXd w(x.size());
    for (int i = 0; i < x.size(); i++) {
        w(i) = x(i).width();
    }
    return w;
}

template <typename Derived>
inline double diagonal_width(const Eigen::MatrixBase<Derived>& x)
{
    Eigen::VectorXd widths = width(x);
    double w = 0;
    for (int i = 0; i < widths.size(); i++) {
        w += widths(i) * widths(i);
    }
    return sqrt(w);
}

template <typename Derived>
inline bool zero_in(const Eigen::MatrixBase<Derived>& x)
{
    // Check if the origin is in the n-dimensional interval
    for (int i = 0; i < x.size(); i++) {
        if (!x(i).zero_in()) {
            return false;
        }
    }
    return true;
}

typedef Vector2<Interval> Vector2I;
typedef Vector3<Interval> Vector3I;
typedef VectorX<Interval> VectorXI;
typedef VectorMax3<Interval> VectorMax3I;
typedef Matrix3<Interval> Matrix2I;
typedef Matrix3<Interval> Matrix3I;
typedef MatrixMax3<Interval> MatrixMax3I;
typedef MatrixX<Interval> MatrixXI;

/// @brief Format a string for an Interval
std::string fmt_interval(const Interval& i, const int precision = 16);
/// @brief Format an eigen VectorX<Interval>
std::string fmt_eigen_intervals(const VectorXI& x, const int precision = 16);
} // namespace ipc::rigid

namespace Eigen {

template <typename BinOp>
struct ScalarBinaryOpTraits<ipc::rigid::Interval, double, BinOp> {
    typedef ipc::rigid::Interval ReturnType;
};

template <typename BinOp>
struct ScalarBinaryOpTraits<double, ipc::rigid::Interval, BinOp> {
    typedef ipc::rigid::Interval ReturnType;
};

// #if EIGEN_MAJOR_VERSION >= 3
// namespace internal {
//     template <typename X, typename S, typename P>
//     struct is_convertible<X, boost::numeric::interval<S, P>> {
//         enum { value = is_convertible<X, S>::value };
//     };

//     template <typename S, typename P1, typename P2>
//     struct is_convertible<
//         boost::numeric::interval<S, P1>,
//         boost::numeric::interval<S, P2>> {
//         enum { value = true };
//     };
// } // namespace internal
// #endif
} // namespace Eigen
