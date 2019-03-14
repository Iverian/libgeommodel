#ifndef GEOM_MODEL_INCLUDE_GM_VEC_H_
#define GEOM_MODEL_INCLUDE_GM_VEC_H_

#include "compare.h"
#include "debug.h"
#include "dot.h"

#include <array>
#include <ostream>

namespace gm {

class Point;

class Vec {
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

    Vec() noexcept;
    Vec(value_type x, value_type y, value_type z) noexcept;
    Vec(value_type magnitude, const std::array<value_type, N>& dir) noexcept;
    Vec(const std::array<value_type, N>& coord) noexcept;
    Vec(const std::initializer_list<value_type>& list) noexcept;

    explicit Vec(const Point& p) noexcept;
    Vec(const Point& lhs, const Point& rhs) noexcept;

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

    Vec& operator+=(const Vec& rhs) noexcept;
    Vec& operator-=(const Vec& rhs) noexcept;
    Vec& operator*=(const_reference rhs) noexcept;
    Vec& operator/=(const_reference rhs) __GM_NOEXCEPT_RELEASE__;

private:
    value_type data_[N];
};

[[nodiscard]] Vec operator-(const Vec& obj) noexcept;

[[nodiscard]] Vec operator+(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] Vec operator-(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] Vec operator*(const Vec& lhs, Vec::const_reference rhs) noexcept;
[[nodiscard]] Vec operator*(Vec::const_reference lhs, const Vec& rhs) noexcept;
[[nodiscard]] Vec operator/(const Vec& lhs,
                            Vec::const_reference rhs) __GM_NOEXCEPT_RELEASE__;

bool operator==(const Vec& lhs, const Vec& rhs) noexcept;
bool operator!=(const Vec& lhs, const Vec& rhs) noexcept;

[[nodiscard]] Vec cross(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] double sqr(const Vec& obj) noexcept;
[[nodiscard]] double norm(const Vec& obj) noexcept;
[[nodiscard]] double dist(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] double cos(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] double sin(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] double angle(const Vec& lhs, const Vec& rhs) noexcept;
[[nodiscard]] Vec unit(const Vec& obj) __GM_NOEXCEPT_RELEASE__;

[[nodiscard]] bool isnan(const Vec& obj) noexcept;
[[nodiscard]] bool isinf(const Vec& obj) noexcept;

std::ostream& operator<<(std::ostream& os, const Vec& obj);

namespace cmp {
    [[nodiscard]] bool near(const Vec& lhs, const Vec& rhs,
                            Tolerance tol = default_tolerance) noexcept;
    [[nodiscard]] bool zero(const Vec& obj,
                            Tolerance tol = default_tolerance) noexcept;
} // namespace cmp

} // namespace gm

namespace std {
template <>
struct hash<gm::Vec> {
    size_t operator()(const gm::Vec& key) const
    {
        return (hasher_(key[0]) << 2) ^ (hasher_(key[1]) << 1)
            ^ hasher_(key[0]);
    }

private:
    hash<gm::Vec::value_type> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_GM_VEC_H_
