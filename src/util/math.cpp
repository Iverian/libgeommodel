#include "math.hpp"

using namespace std;

double pad(double t, double f, double b, double pad)
{
    return (b * (t + pad) - f * (t - pad) - 2 * pad * t) / (b - f);
}

array<double, 2> atan2v(double y, double x)
{
    auto u1 = atan2(y, x);
    if (u1 < 0) {
        u1 += M_PI;
    }
    auto u2 = u1 + M_PI;
    return {u1, u2};
}

double binom(size_t n, size_t k) __GM_NOEXCEPT_RELEASE__
{
    check_ifd(k <= n, "Unable to compute C^n_k with n = {}, k = {}", n, k);

    double result = 1;
    if (k != 0 && k != n) {
        for (size_t i = 0; i < k; ++i) {
            result *= double(n - i) / (k - i);
        }
    }
    return result;
}