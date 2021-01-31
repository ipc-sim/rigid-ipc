#pragma once

#if 1
namespace ipc::rigid {

class Multiprecision {
private:
    double _var;

public:
    Multiprecision(const double var = 0, int dummy = 0)
        : _var(var)
    {
    }

    operator double() { return _var; }

    Multiprecision& operator=(const double& var)
    {
        _var = var;
        return *this;
    }
};
} // namespace ipc::rigid
#else
#include <cmath>
#include <gmp.h>
#include <iostream>

namespace ipc::rigid {

class Multiprecision {
public:
    static void set_precision(const int digits)
    {
        mpf_set_default_prec(digits);
    }

    static int get_precision() { return mpf_get_default_prec(); }

    mpf_t value;
    bool isinf = false;

    int get_sign() { return mpf_sgn(value); }
    int get_prec_bits() { return mpf_get_prec(value); }
    Multiprecision()
    {
        mpf_init(value);
        mpf_set_d(value, 0);
    }

    Multiprecision(double d)
    {
        isinf = std::isinf(d);
        mpf_init(value);
        mpf_set_d(value, isinf ? 10E100 : d);
    }

    Multiprecision(double d, int precision)
    {
        isinf = std::isinf(d);
        mpf_init2(value, precision);
        mpf_set_d(value, isinf ? 10E100 : d);
    }

    Multiprecision(const mpf_t& v_)
    {
        mpf_init2(value, mpf_get_prec(v_));
        mpf_set(value, v_);
    }

    Multiprecision(const Multiprecision& other)
    {
        isinf = other.isinf;
        mpf_init2(value, mpf_get_prec(other.value));
        mpf_set(value, other.value);
    }

    ~Multiprecision() { mpf_clear(value); }

    friend Multiprecision
    operator+(const Multiprecision& x, const Multiprecision& y)
    {
        static Multiprecision r_out;
        mpf_add(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Multiprecision
    operator-(const Multiprecision& x, const Multiprecision& y)
    {
        static Multiprecision r_out;
        mpf_sub(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Multiprecision
    operator*(const Multiprecision& x, const Multiprecision& y)
    {
        static Multiprecision r_out;
        mpf_mul(r_out.value, x.value, y.value);
        return r_out;
    }

    friend Multiprecision
    operator/(const Multiprecision& x, const Multiprecision& y)
    {
        static Multiprecision r_out;
        mpf_div(r_out.value, x.value, y.value);
        return r_out;
    }

    Multiprecision& operator=(const Multiprecision& x)
    {
        if (this == &x)
            return *this;
        // mpf_init2(value, prec);
        mpf_set(value, x.value);
        isinf = x.isinf;
        return *this;
    }

    Multiprecision& operator=(const double x)
    {
        // mpf_init2(value, prec);
        mpf_set_d(value, x);
        return *this;
    }

    //> < ==
    friend bool operator<(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) < 0;
    }

    friend bool operator>(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) > 0;
    }

    friend bool operator<=(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) <= 0;
    }

    friend bool operator>=(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) >= 0;
    }

    friend bool operator==(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) == 0;
    }

    friend bool operator!=(const Multiprecision& r, const Multiprecision& r1)
    {
        return mpf_cmp(r.value, r1.value) != 0;
    }

    // to double
    double to_double()
    {
        return (isinf || mpf_get_d(value) > 10E99)
            ? std::numeric_limits<double>::infinity()
            : mpf_get_d(value);
    }

    //<<
    friend std::ostream& operator<<(std::ostream& os, const Multiprecision& r)
    {
        os << mpf_get_d(r.value);
        return os;
    }

    friend Multiprecision sqrt(const Multiprecision& mp)
    {
        Multiprecision res;
        mpf_sqrt(res.value, mp.value);

        return res;
    }
};
} // namespace ipc::rigid
#endif
