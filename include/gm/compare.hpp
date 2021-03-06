#ifndef GEOM_MODEL_INCLUDE_GM_COMPARE_HPP__
#define GEOM_MODEL_INCLUDE_GM_COMPARE_HPP__

#include "exports.hpp"

#include <ostream>

namespace gm {

enum class Tolerance : unsigned {
    ZERO = 0,
    SINGLE = 1,
    DOUBLE = 2,
    TRIPLE = 3,
    MAX = 4
};

GM_EXPORT std::ostream& operator<<(std::ostream& os, const Tolerance& obj);

namespace cmp {
    static constexpr auto default_tolerance = Tolerance::DOUBLE;

    [[nodiscard]] GM_EXPORT double tol(Tolerance t
                                       = default_tolerance) noexcept;

    [[nodiscard]] GM_EXPORT bool
    isnear(double lhs, double rhs, Tolerance t = default_tolerance) noexcept;

    [[nodiscard]] GM_EXPORT bool
    near(double lhs, double rhs, Tolerance t = default_tolerance) noexcept;

    [[nodiscard]] GM_EXPORT bool
    zero(double lhs, Tolerance t = default_tolerance) noexcept;

    [[nodiscard]] GM_EXPORT bool le(double lhs, double rhs,
                                    Tolerance t = default_tolerance) noexcept;
    [[nodiscard]] GM_EXPORT bool ge(double lhs, double rhs,
                                    Tolerance t = default_tolerance) noexcept;
} // namespace cmp

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_COMPARE_HPP__
