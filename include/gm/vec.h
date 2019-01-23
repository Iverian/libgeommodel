#ifndef GEOM_MODEL_INCLUDE_VEC_H_
#define GEOM_MODEL_INCLUDE_VEC_H_

#include <array>
#include <ostream>

#include "compare.h"

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

    pointer data() noexcept;
    const_pointer data() const noexcept;
    size_type size() const noexcept;
    reference operator[](size_type i) noexcept;
    const_reference operator[](size_type i) const noexcept;

    iterator begin() noexcept;
    iterator end() noexcept;
    const_iterator begin() const noexcept;
    const_iterator end() const noexcept;

    std::array<value_type, N> raw() const noexcept;

    Vec& operator+=(const Vec& rhs) noexcept;
    Vec& operator-=(const Vec& rhs) noexcept;
    Vec& operator*=(const_reference rhs) noexcept;
    Vec& operator/=(const_reference rhs) noexcept;

private:
    value_type data_[N];
};

Vec operator-(const Vec& obj) noexcept;

Vec operator+(const Vec& lhs, const Vec& rhs) noexcept;
Vec operator-(const Vec& lhs, const Vec& rhs) noexcept;
Vec operator*(const Vec& lhs, Vec::const_reference rhs) noexcept;
Vec operator*(Vec::const_reference lhs, const Vec& rhs) noexcept;
Vec operator/(const Vec& lhs, Vec::const_reference rhs) noexcept;

bool operator==(const Vec& lhs, const Vec& rhs) noexcept;
bool operator!=(const Vec& lhs, const Vec& rhs) noexcept;

double dot(const Vec& lhs, const Vec& rhs) noexcept;
Vec cross(const Vec& lhs, const Vec& rhs) noexcept;
double sqr(const Vec& obj) noexcept;
double norm(const Vec& obj) noexcept;
double dist(const Vec& lhs, const Vec& rhs) noexcept;
double cos(const Vec& lhs, const Vec& rhs) noexcept;
double sin(const Vec& lhs, const Vec& rhs) noexcept;
double angle(const Vec& lhs, const Vec& rhs) noexcept;
Vec unit(const Vec& obj) noexcept;

bool isnan(const Vec& obj) noexcept;
bool isinf(const Vec& obj) noexcept;
bool isnear(const Vec& lhs, const Vec& rhs,
            Tolerance tol = Tolerance::DOUBLE) noexcept;

std::ostream& operator<<(std::ostream& os, const Vec& obj);

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

#endif // GEOM_MODEL_INCLUDE_VEC_H_
