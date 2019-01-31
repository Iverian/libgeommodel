#ifndef GEOM_MODEL_INCLUDE_DOT_H_
#define GEOM_MODEL_INCLUDE_DOT_H_

#include <iterator>
#include <type_traits>

namespace gm {

template <class T, class U,
          typename = std::enable_if_t<
              std::is_same_v<typename T::value_type, typename U::value_type>>>
[[nodiscard]] typename T::value_type dot(const T& lhs, const U& rhs) noexcept
{
    using value_type = typename T::value_type;

    value_type result {};
    auto i = std::begin(lhs);
    auto ilast = std::end(lhs);
    auto j = std::begin(rhs);
    auto jlast = std::end(rhs);

    for (;i != ilast && j != jlast; ++i, ++j) {
        result += (*i) * (*j);
    }
    return result;
}

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_DOT_H_