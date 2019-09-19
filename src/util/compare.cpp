#include <cmath>
#include <gm/compare.hpp>

namespace gm {

namespace cmp {

double tol(Tolerance t) noexcept
{
    static constexpr double tol[] = {1e-2, 1e-5, 1e-10, 1e-25, 1e-35};
    return tol[size_t(t)];
}

bool isnear(double lhs, double rhs, Tolerance t) noexcept
{
    return fabs(rhs - lhs) < tol(t);
}

bool near(double lhs, double rhs, Tolerance t) noexcept
{
    return isnear(lhs, rhs, t);
}

bool zero(double lhs, Tolerance t) noexcept
{
    return fabs(lhs) < tol(t);
}

bool le(double lhs, double rhs, Tolerance t) noexcept
{
    return lhs < rhs || isnear(lhs, rhs, t);
}

bool ge(double lhs, double rhs, Tolerance t) noexcept
{
    return lhs > rhs || isnear(lhs, rhs, t);
}

} // namespace cmp

std::ostream& operator<<(std::ostream& os, const Tolerance& obj)
{
    return os << cmp::tol(obj);
}

} // namespace gm
