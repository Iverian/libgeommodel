#ifndef GEOM_MODEL_INCLUDE_GM_POINT_HPP_
#define GEOM_MODEL_INCLUDE_GM_POINT_HPP_

#include "compare.hpp"
#include "debug.hpp"
#include "dot.hpp"
#include "exports.hpp"

#include <array>
#include <ostream>

namespace gm {

class GM_EXPORT Vec;

class GM_EXPORT Point {
public:
    static constexpr size_t N = 3;
    using value_type = double;
    using size_type = size_t;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using iterator = pointer;
    using const_iterator = const_pointer;

    Point() noexcept;
    Point(value_type x, value_type y, value_type z) noexcept;
    Point(const std::array<value_type, N>& coord) noexcept;
    Point(const std::initializer_list<value_type>& list) noexcept;

    Point(const Vec& lhs, const Vec& rhs) noexcept;
    explicit Point(const Vec& v) noexcept;

    [[nodiscard]] pointer data() noexcept;
    [[nodiscard]] const_pointer data() const noexcept;
    [[nodiscard]] size_type size() const noexcept;
    [[nodiscard]] reference operator[](size_type i) noexcept;
    [[nodiscard]] const_reference operator[](size_type i) const noexcept;

    [[nodiscard]] iterator begin() noexcept;
    [[nodiscard]] iterator end() noexcept;
    [[nodiscard]] const_iterator begin() const noexcept;
    [[nodiscard]] const_iterator end() const noexcept;

    [[nodiscard]] std::array<value_type, N> raw() const noexcept;

    Point& operator+=(const Point& rhs) noexcept;
    Point& operator-=(const Point& rhs) noexcept;
    Point& operator*=(const_reference rhs) noexcept;
    Point& operator/=(const_reference rhs) __GM_NOEXCEPT_RELEASE__;

private:
    value_type data_[N];
};

[[nodiscard]] GM_EXPORT Point operator-(const Point& obj) noexcept;

[[nodiscard]] GM_EXPORT Point operator+(const Point& lhs,
                                        const Point& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator-(const Point& lhs,
                                        const Point& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator*(const Point& lhs,
                                        Point::const_reference rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator*(Point::const_reference lhs,
                                        const Point& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator/(
    const Point& lhs, Point::const_reference rhs) __GM_NOEXCEPT_RELEASE__;

[[nodiscard]] GM_EXPORT Point operator+(const Point& lhs,
                                        const Vec& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator+(const Vec& lhs,
                                        const Point& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator-(const Point& lhs,
                                        const Vec& rhs) noexcept;
[[nodiscard]] GM_EXPORT Point operator-(const Vec& lhs,
                                        const Point& rhs) noexcept;

[[nodiscard]] GM_EXPORT bool operator==(const Point& lhs,
                                        const Point& rhs) noexcept;
[[nodiscard]] GM_EXPORT bool operator!=(const Point& lhs,
                                        const Point& rhs) noexcept;

[[nodiscard]] GM_EXPORT double sqr(const Point& obj) noexcept;
[[nodiscard]] GM_EXPORT double dist(const Point& lhs,
                                    const Point& rhs) noexcept;

[[nodiscard]] GM_EXPORT bool isnan(const Point& obj) noexcept;
[[nodiscard]] GM_EXPORT bool isinf(const Point& obj) noexcept;

GM_EXPORT std::ostream& operator<<(std::ostream& os, const Point& obj);

namespace cmp {
    [[nodiscard]] GM_EXPORT bool near(const Point& lhs, const Point& rhs,
                                      Tolerance tol
                                      = default_tolerance) noexcept;
    [[nodiscard]] GM_EXPORT bool
    zero(const Point& obj, Tolerance tol = default_tolerance) noexcept;
} // namespace cmp

} // namespace gm

namespace std {
template <>
struct hash<gm::Point> {
    size_t operator()(const gm::Point& key) const
    {
        return (hasher_(key[0]) << 2) ^ (hasher_(key[1]) << 1)
            ^ hasher_(key[2]);
    }

private:
    hash<gm::Point::value_type> hasher_;
};
} // namespace std

#endif // GEOM_MODEL_INCLUDE_GM_POINT_HPP_
