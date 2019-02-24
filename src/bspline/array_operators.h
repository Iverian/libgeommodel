#ifndef GEOM_MODEL_SRC_BSPLINE_ARRAY_OPERATORS_H_
#define GEOM_MODEL_SRC_BSPLINE_ARRAY_OPERATORS_H_

#include <array>
#include <type_traits>

#include <gm/compare.h>
#include <util/debug.h>

namespace gm {

template <class T, class U, size_t N>
std::array<T, N>& operator+=(std::array<T, N>& lhs,
                             const std::array<U, N>& rhs) noexcept
{
    for (size_t i = 0; i < N; ++i) {
        lhs[i] += rhs[i];
    }
    return lhs;
}

template <class T, class U, size_t N>
std::array<T, N>& operator-=(std::array<T, N>& lhs,
                             const std::array<U, N>& rhs) noexcept
{
    for (size_t i = 0; i < N; ++i) {
        lhs[i] -= rhs[i];
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N>&
operator*=(std::array<T, N>& lhs,
           const typename std::array<T, N>::const_reference rhs) noexcept
{
    for (size_t i = 0; i < N; ++i) {
        lhs[i] *= rhs;
    }
    return lhs;
}

template <class T, size_t N>
std::array<T, N>&
operator/=(std::array<T, N>& lhs,
           const typename std::array<T, N>::const_reference rhs)
    __GM_NOEXCEPT_RELEASE__
{
    check_ifd(!gm::cmp::zero(rhs), "Division by zero");

    for (size_t i = 0; i < N; ++i) {
        lhs[i] /= rhs;
    }
    return lhs;
}

template <class T, class U, size_t N>
std::array<std::common_type_t<T, U>, N>
operator+(const std::array<T, N>& lhs, const std::array<U, N>& rhs) noexcept
{
    std::array<std::common_type_t<T, U>, N> result {};
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

template <class T, class U, size_t N>
std::array<std::common_type_t<T, U>, N>
operator-(const std::array<T, N>& lhs, const std::array<U, N>& rhs) noexcept
{
    std::array<std::common_type_t<T, U>, N> result {};
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

template <class T, size_t N>
std::array<T, N>
operator*(const std::array<T, N>& lhs,
          const typename std::array<T, N>::const_reference rhs) noexcept
{
    std::array<T, N> result {};
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] * rhs;
    }
    return result;
}

template <class T, size_t N>
std::array<T, N>
operator*(const typename std::array<T, N>::const_reference lhs,
          const std::array<T, N>& rhs) noexcept
{
    std::array<T, N> result {};
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs * rhs[i];
    }
    return result;
}

template <class T, size_t N>
std::array<T, N>
operator/(const std::array<T, N>& lhs,
          const typename std::array<T, N>::const_reference rhs)
    __GM_NOEXCEPT_RELEASE__
{
    check_ifd(!gm::cmp::zero(rhs), "Division by zero");

    std::array<T, N> result {};
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] / rhs;
    }
    return result;
}

} // namespace gm

#endif // GEOM_MODEL_SRC_BSPLINE_ARRAY_OPERATORS_H_