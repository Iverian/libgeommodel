#include <cmath>
#include <gm/compare.h>

using namespace std;

namespace gm {

double tol(Tolerance t) noexcept
{
    static constexpr double tol[] = {1e-2, 1e-5, 1e-9, 1e-20};
    return tol[size_t(t)];
}

bool isnear(double lhs, double rhs, Tolerance t) noexcept
{

    return fabs(rhs - lhs) < tol(t);
}

bool iszero(double lhs, Tolerance t) noexcept
{
    return isnear(lhs, 0., t);
}

} // namespace gm