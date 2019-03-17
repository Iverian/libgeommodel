#include <cmath>
#include <gm/compare.hpp>

namespace gm::cmp {
double tol(Tolerance t) noexcept
{
    static constexpr double tol[] = {1e-2, 1e-5, 1e-9, 1e-20};
    return tol[size_t(t)];
}

bool near(double lhs, double rhs, Tolerance t) noexcept
{

    return fabs(rhs - lhs) < tol(t);
}

bool zero(double lhs, Tolerance t) noexcept
{
    return fabs(lhs) < tol(t);
}

bool le(double lhs, double rhs, Tolerance t) noexcept
{
    return lhs < rhs || near(lhs, rhs, t);
}

bool ge(double lhs, double rhs, Tolerance t) noexcept
{
    return lhs > rhs || near(lhs, rhs, t);
}

} // namespace gm::cmp