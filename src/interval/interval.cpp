#include "interval.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>

namespace ipc::rigid {

Interval::Interval() { *this = empty(); }

Interval::Interval(double x)
{
    INF = x;
    SUP = x;
}

Interval::Interval(double x, double y)
{
    INF = x;
    SUP = y;
}

Interval Interval::empty()
{
    return Interval(
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity());
}

/* ------------------------------------------------------------------- */
/* --- IO (input/output)                                           --- */
/* ------------------------------------------------------------------- */

std::istream& operator>>(std::istream& is, Interval& a)
{
    double help, ioconst;

    ioconst = (1e17 - 1) * 1e27;

    is >> help;
    if ((help < -ioconst) || (help > ioconst))
        a.INF = q_pred(q_pred(help));
    else
        a.INF = q_pred(help);

    is >> help;
    if ((help < -ioconst) || (help > ioconst))
        a.SUP = q_succ(q_succ(help));
    else
        a.SUP = q_succ(help);

    return is;
}

std::ostream& operator<<(std::ostream& os, const Interval& a)
{
    Interval help;

    help.INF = q_pred(q_pred(a.INF));
    help.SUP = q_succ(q_succ(a.SUP));

    std::ios_base::fmtflags aktform = std::cout.flags();
    os << "[" << std::setprecision(15) << std::setw(23)
       << std::setiosflags(std::ios::scientific);
    os << help.INF;
    std::cout.flags(aktform);
    os << "," << std::setprecision(15) << std::setw(23)
       << std::setiosflags(std::ios::scientific);
    os << help.SUP;
    std::cout.flags(aktform);
    os << " ]";

    return os;
}

/* ------------------------------------------------------------------- */
/* --- interval arithmetic (basic operations)                      --- */
/* ------------------------------------------------------------------- */

Interval operator+(const Interval& a, const Interval& b)
{
    return add_ii(a, b);
}

Interval operator+(const Interval& a, double b) { return add_id(a, b); }

Interval operator+(double a, const Interval& b) { return add_di(a, b); }

Interval operator+(const Interval& a) { return a; }

Interval& Interval::operator+=(const Interval& rhs)
{
    *this = add_ii(*this, rhs);
    return *this;
}

Interval& Interval::operator+=(const double& rhs)
{
    *this = add_id(*this, rhs);
    return *this;
}

Interval operator-(const Interval& a, const Interval& b)
{
    return sub_ii(a, b);
}

Interval operator-(const Interval& a, double b) { return sub_id(a, b); }

Interval operator-(double a, const Interval& b) { return sub_di(a, b); }

Interval operator-(const Interval& a) { return Interval(-a.SUP, -a.INF); }

Interval& Interval::operator-=(const Interval& rhs)
{
    *this = sub_ii(*this, rhs);
    return *this;
}

Interval& Interval::operator-=(const double rhs)
{
    *this = sub_id(*this, rhs);
    return *this;
}

Interval operator*(const Interval& a, const Interval& b)
{
    return mul_ii(a, b);
}

Interval operator*(const Interval& a, double b) { return mul_id(a, b); }

Interval operator*(double a, const Interval& b) { return mul_di(a, b); }

Interval& Interval::operator*=(const Interval& rhs)
{
    *this = mul_ii(*this, rhs);
    return *this;
}

Interval& Interval::operator*=(const double rhs)
{
    *this = mul_id(*this, rhs);
    return *this;
}

Interval operator/(const Interval& a, const Interval& b)
{
    return div_ii(a, b);
}

Interval operator/(const Interval& a, double b) { return div_id(a, b); }

Interval operator/(double a, const Interval& b) { return div_di(a, b); }

Interval& Interval::operator/=(const Interval& rhs)
{
    *this = div_ii(*this, rhs);
    return *this;
}

Interval& Interval::operator/=(const double& rhs)
{
    *this = div_id(*this, rhs);
    return *this;
}

/* ------------------------------------------------------------------- */
/* --- Interval arithmetic (logical operations)                    --- */
/* ------------------------------------------------------------------- */

bool operator==(const Interval& a, const Interval& b) { return ieq_ii(a, b); }

bool operator==(const Interval& a, double b) { return ieq_ii(a, Interval(b)); }

bool operator!=(const Interval& a, const Interval& b)
{
    return !(a.INF == b.INF && a.SUP == b.SUP);
}

bool operator<(const Interval& a, const Interval& b) { return in_ii(a, b); }

bool operator<(const double a, const Interval& b)
{
    return b.INF < a && a < b.SUP;
}

bool operator<=(const Interval& a, const Interval& b)
{
    return b.INF <= a.INF && a.SUP <= b.SUP;
}

bool operator<=(double a, const Interval& b) { return in_di(a, b); }

bool operator>(const Interval& a, const Interval& b)
{
    return b.INF > a.INF && a.SUP > b.SUP;
}
bool operator>(const Interval& a, double b) { return a.INF < b && b < a.SUP; }

bool operator>=(const Interval& a, const Interval& b)
{
    return b.INF >= a.INF && a.SUP >= b.SUP;
}

bool operator>=(const Interval& a, double b)
{
    return a.INF <= b && b <= a.SUP;
}

Interval operator|(const Interval& a, const Interval& b) { return hull(a, b); }
Interval& Interval::operator|=(const Interval& b)
{
    *this = *this | b;
    return *this;
}

Interval operator&(const Interval& a, const Interval& b)
{
    return intsec(a, b);
}
Interval& Interval::operator&=(const Interval& b)
{
    *this = *this & b;
    return *this;
}

bool in(double a, const Interval& b) { return in_di(a, b); }

bool in(const Interval& a, const Interval& b) { return in_ii(a, b); }

bool Interval::is_empty() const { return this->INF > this->SUP; }

bool Interval::zero_in() const { return this->INF <= 0 && this->SUP >= 0; }

bool Interval::overlaps(const Interval& b) const
{
    return (this->INF <= b.SUP) && (b.INF <= this->SUP);
}

/* ------------------------------------------------------------------- */
/* --- utilities, mid, diam, ...                                   --- */
/* ------------------------------------------------------------------- */

double Interval::mid() const { return q_mid(*this); }

std::pair<Interval, Interval> Interval::bisect() const
{
    double mid = this->mid();
    return std::make_pair(Interval(this->INF, mid), Interval(mid, this->SUP));
}

bool disjoint(const Interval& a, const Interval& b) { return dis_ii(a, b); }

double Interval::diam() const { return q_diam(*this); }

double Interval::drel() const
{
    if ((SUP <= -q_minr) || (q_minr <= INF)) {
        if (INF > 0)
            return diam() / INF;
        else
            return diam() / (-SUP);
    } else {
        return diam();
    }
}

Interval Interval::blow(double eps) const
{
    Interval y;
    y = (1.0 + eps) * (*this) - eps * (*this);
    return (Interval(q_pred(y.INF), q_succ(y.SUP)));
}

// min max function, same as what BOOST does
Interval max(const Interval& x, const Interval& y)
{
    return Interval(std::max(x.INF, y.INF), std::max(x.SUP, y.SUP));
}
Interval max(const Interval& x, double y)
{
    return Interval(std::max(x.INF, y), std::max(x.SUP, y));
}

Interval max(double x, const Interval& y)
{
    return Interval(std::max(x, y.INF), std::max(x, y.SUP));
}

Interval min(const Interval& x, const Interval& y)
{
    return Interval(std::min(x.INF, y.INF), std::min(x.SUP, y.SUP));
}
Interval min(const Interval& x, double y)
{
    return Interval(std::min(x.INF, y), std::min(x.SUP, y));
}

Interval min(double x, const Interval& y)
{
    return Interval(std::min(x, y.INF), std::min(x, y.SUP));
}

/* ------------------------------------------------------------------- */
/* --- Interval arithmetic (elementary functions)                  --- */
/* ------------------------------------------------------------------- */

Interval exp(const Interval& a) { return j_exp(a); }
Interval expm(const Interval& a) { return j_expm(a); }
Interval sinh(const Interval& a) { return j_sinh(a); }
Interval cosh(const Interval& a) { return j_cosh(a); }
Interval coth(const Interval& a) { return j_coth(a); }
Interval tanh(const Interval& a) { return j_tanh(a); }
Interval log(const Interval& a) { return j_log(a); }
Interval ln(const Interval& a) { return j_log(a); }
Interval lg1p(const Interval& a) { return j_lg1p(a); }
Interval sqrt(const Interval& a) { return j_sqrt(a); }
Interval sqr(const Interval& a) { return j_sqr(a); }
Interval asnh(const Interval& a) { return j_asnh(a); }
Interval asinh(const Interval& a) { return j_asnh(a); }
Interval acsh(const Interval& a) { return j_acsh(a); }
Interval acosh(const Interval& a) { return j_acsh(a); }
Interval acth(const Interval& a) { return j_acth(a); }
Interval acoth(const Interval& a) { return j_acth(a); }
Interval atnh(const Interval& a) { return j_atnh(a); }
Interval atanh(const Interval& a) { return j_atnh(a); }
Interval asin(const Interval& a) { return j_asin(a); }
Interval acos(const Interval& a) { return j_acos(a); }
Interval acot(const Interval& a) { return j_acot(a); }
Interval atan(const Interval& a) { return j_atan(a); }
Interval sin(const Interval& a) { return j_sin(a); }
Interval cos(const Interval& a) { return j_cos(a); }
Interval cot(const Interval& a) { return j_cot(a); }
Interval tan(const Interval& a) { return j_tan(a); }
Interval exp2(const Interval& a) { return j_exp2(a); }
Interval ex10(const Interval& a) { return j_ex10(a); }
Interval log2(const Interval& a) { return j_log2(a); }
Interval lg10(const Interval& a) { return j_lg10(a); }
Interval erf(const Interval& a) { return j_erf(a); }
Interval erfc(const Interval& a) { return j_erfc(a); }

Interval abs(const Interval& a)
{
    if (a.INF >= 0)
        return a;
    else if (a.SUP <= 0)
        return -a;
    else
        return Interval(0, std::max(-a.INF, a.SUP));
}

} // namespace ipc::rigid