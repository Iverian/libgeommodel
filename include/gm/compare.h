#ifndef GEOM_MODEL_INCLUDE_GM_COMPARE_H__
#define GEOM_MODEL_INCLUDE_GM_COMPARE_H__

namespace gm {

enum class Tolerance : unsigned { ZERO = 0, SINGLE = 1, DOUBLE = 2, MAX = 3 };

static constexpr auto default_tolerance = Tolerance::DOUBLE;

[[nodiscard]] double tol(Tolerance t) noexcept;
[[nodiscard]] bool isnear(double lhs, double rhs,
                          Tolerance t = default_tolerance) noexcept;
[[nodiscard]] bool iszero(double lhs,
                          Tolerance t = default_tolerance) noexcept;

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_COMPARE_H__