#ifndef GEOM_MODEL_INCLUDE_DOT_H_
#define GEOM_MODEL_INCLUDE_DOT_H_

#include <iterator>
#include <type_traits>

namespace gm {

template <class T, class U>
using dottable_t = std::enable_if_t<
    std::is_arithmetic_v<
        T> && std::is_arithmetic_v<U> && std::is_convertible_v<T, U> && std::is_convertible_v<U, T>>;

template <class T, class U,
          class = dottable_t<typename T::value_type, typename U::value_type>>
[[nodiscard]] std::common_type_t<typename T::value_type,
                                 typename U::value_type>
dot(const T& lhs, const U& rhs) noexcept
{
    std::common_type_t<typename T::value_type, typename U::value_type>
        result {};
    auto i = std::begin(lhs);
    auto ilast = std::end(lhs);
    auto j = std::begin(rhs);
    auto jlast = std::end(rhs);

    for (; i != ilast && j != jlast; ++i, ++j) {
        result += (*i) * (*j);
    }
    return result;
}
template <class T,
          class = dottable_t<typename T::value_type, typename T::value_type>>
[[nodiscard]] typename T::value_type sqr(const T& t) noexcept
{
    return dot(t, t);
}

template <class T, class = std::enable_if_t<std::is_arithmetic_v<T>>>
[[nodiscard]] T sqr(const T& t)
{
    return t * t;
}

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_DOT_H_