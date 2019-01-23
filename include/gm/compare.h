#ifndef GEOM_MODEL_INCLUDE_GEOM_MODEL_COMPARE_H_
#define GEOM_MODEL_INCLUDE_GEOM_MODEL_COMPARE_H_

namespace gm {

enum class Tolerance : unsigned { ZERO = 0, SINGLE = 1, DOUBLE = 2, MAX = 3 };

double tol(Tolerance t) noexcept;
bool isnear(double lhs, double rhs, Tolerance t = Tolerance::DOUBLE) noexcept;
bool iszero(double lhs, Tolerance t = Tolerance::DOUBLE) noexcept;

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GEOM_MODEL_COMPARE_H_