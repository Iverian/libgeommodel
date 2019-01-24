#ifndef GEOM_MODEL_INCLUDE_GM_POINT_H_
#define GEOM_MODEL_INCLUDE_GM_POINT_H_

#include <array>
#include <ostream>

#include "compare.h"

namespace gm {

class Vec;

class Point {
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

    Point& operator+=(const Point& rhs) noexcept;
    Point& operator-=(const Point& rhs) noexcept;
    Point& operator*=(const_reference rhs) noexcept;
    Point& operator/=(const_reference rhs) noexcept;

private:
    value_type data_[N];
};

Point operator-(const Point& obj) noexcept;

Point operator+(const Point& lhs, const Point& rhs) noexcept;
Point operator-(const Point& lhs, const Point& rhs) noexcept;
Point operator*(const Point& lhs, Point::const_reference rhs) noexcept;
Point operator*(Point::const_reference lhs, const Point& rhs) noexcept;
Point operator/(const Point& lhs, Point::const_reference rhs) noexcept;

Point operator+(const Point& lhs, const Vec& rhs) noexcept;
Point operator+(const Vec& lhs, const Point& rhs) noexcept;
Point operator-(const Point& lhs, const Vec& rhs) noexcept;
Point operator-(const Vec& lhs, const Point& rhs) noexcept;

bool operator==(const Point& lhs, const Point& rhs) noexcept;
bool operator!=(const Point& lhs, const Point& rhs) noexcept;

double dot(const Point& lhs, const Point& rhs) noexcept;
double dot(const Vec& lhs, const Point& rhs) noexcept;
double dot(const Point& lhs, const Vec& rhs) noexcept;

double sqr(const Point& obj) noexcept;
double dist(const Point& lhs, const Point& rhs) noexcept;

bool isnan(const Point& obj) noexcept;
bool isinf(const Point& obj) noexcept;
bool isnear(const Point& lhs, const Point& rhs,
            Tolerance tol = Tolerance::DOUBLE) noexcept;

std::ostream& operator<<(std::ostream& os, const Point& obj);

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
}

#endif // GEOM_MODEL_INCLUDE_GM_POINT_H_
