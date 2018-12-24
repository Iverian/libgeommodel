#ifndef GEOM_MODEL_INCLUDE_VEC_H_
#define GEOM_MODEL_INCLUDE_VEC_H_

#include "dot.h"

#include <array>
#include <memory>

class Vec;
class Point;
class VecImpl;

double norm(const Vec& x);
double dist(const Vec& x, const Vec& y);
double sqr(const Vec& a);
double cos(const Vec& a, const Vec& b);
double sin(const Vec& a, const Vec& b);
double angle(const Vec& a, const Vec& b);
Vec cross(const Vec& a, const Vec& b);
Vec operator*(double x, const Vec& v);

Vec unit(const Vec& v);
Vec unit(Vec&& v);

class Vec {
    friend class Point;

public:
    Vec();
    Vec(double x, double y, double z);
    Vec(double magnitude, const std::array<double, 3>& coord);
    Vec(const std::initializer_list<double>& list);
    Vec(const Point& a, const Point& b);
    Vec(const std::pair<Point, Point>& diff);
    explicit Vec(const std::array<double, 3>& coord);
    explicit Vec(const Point& p);
    explicit Vec(Point&& p);

    std::array<double, 3>::iterator begin();
    std::array<double, 3>::iterator end();
    std::array<double, 3>::const_iterator cbegin() const;
    std::array<double, 3>::const_iterator cend() const;

    const double& operator[](size_t i) const noexcept;
    size_t size() const noexcept;
    const std::array<double, 3>& raw() const noexcept;

    Vec& operator+=(const Vec& other);
    Vec& operator-=(const Vec& other);
    Vec& operator*=(double x);
    Vec& operator/=(double x);
    Vec operator+(const Vec& other) const;
    Vec operator-() const;
    Vec operator-(const Vec& other) const;
    Vec operator*(double x) const;
    Vec operator/(double x) const;

    double sqr() const;
    double norm() const;
    bool isnan() const;

    Vec& normalize();
    Vec normalize() const;

    friend bool operator==(const Vec& lhs, const Vec& rhs);
    friend bool operator!=(const Vec& lhs, const Vec& rhs);
    friend std::ostream& operator<<(std::ostream& os, const Vec& v);

private:
    void copy();

    std::shared_ptr<VecImpl> pimpl_;
};

namespace std {
template <class T, size_t N>
struct hash<array<T, N>> {
    static constexpr size_t seed = 13873;
    using argument_type = array<T, N>;

    size_t operator()(const argument_type& key) const
    {
        size_t result = seed;
        for (size_t i = 0; i < key.size(); ++i) {
            result = (result << 1) ^ hasher_(key[i]);
        }
        return result;
    }

private:
    hash<T> hasher_;
};

template <>
struct hash<Vec> {
    size_t operator()(const Vec& key) const
    {
        return hasher_(key.raw());
    }

private:
    hash<array<double, 3>> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_VEC_H_
