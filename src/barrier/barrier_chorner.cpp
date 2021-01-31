#include <iomanip>
#include <iostream>

#include <cmath>
#include <limits>

#include <assert.h>
#include <math.h>

namespace ipc::rigid {

double horner_double(double* polynomial, int n, double x)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s = s * x + polynomial[i];
    }
    return s;
}

double horner_fma(double* polynomial, int n, double x)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) {
        s = fma(s, x, polynomial[i]);
    }
    return s;
}

void _two_sum(double a, double b, double& x, double& y)
{
    double t;
    x = a + b;
    t = x - a;
    y = (a - (x - t)) + (b - t);
}

void _two_product_fma(double a, double b, double& x, double& y)
{
    x = a * b;
    y = fma(a, b, -x);
}

double horner_compensated(double* polynomial, int n, double x)
{
    // Based on "Compensated Horner Scheme"
    // by S. Graillat, Ph. Langlois, N. Louvet, 2005
    double s = 0.0;
    double err = 0.0;
    for (int i = 0; i < n; i++) {
        double p, pi, sigma;
        _two_product_fma(s, x, p, pi);
        _two_sum(p, polynomial[i], s, sigma);
        err = err * x + (pi + sigma);
    }
    return s + err;
}

double barrier_horner_compensated(double x, double eps)
{
    assert(eps > 0);

    if (x <= 0.0)
        return std::numeric_limits<double>::infinity();

    if (x >= eps)
        return 0.0;

    double y = x / eps;
    double numer = (1.0 - y);
    numer = numer * numer * numer;
    double denom_coeff[] = { 1.0, -3.0, 3.0, 0.0 };
    double denom = horner_compensated(denom_coeff, 4, y);

    return numer / denom;
}

} // namespace ipc::rigid
