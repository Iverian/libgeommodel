#ifndef GEOM_MODEL_INCLUDE_POINT_H_
#define GEOM_MODEL_INCLUDE_POINT_H_

#include "vec.h"

#include <memory>

class Point;
class VecImpl;

double dist(const Point& a, const Point& b);
bool isnear(const Point& lhs, const Point& rhs, double eps = 1e-5);
double angle(const Point& a, const Point& b, const Point& c, const Vec& norm);
Point operator+(const Point& p, const Vec& v);
Point operator+(const Vec& v, const Point& p);
Point operator-(const Point& p, const Vec& v);
Point operator-(const Vec& v, const Point& p);

Point operator*(double x, const Point& p);

class Point {
    friend class Vec;

public:
    Point();
    Point(double x, double y, double z);
    Point(const std::initializer_list<double>& list);
    explicit Point(const std::array<double, 3>& coord);
    explicit Point(const Vec& v);
    explicit Point(Vec&& v);

    std::array<double, 3>::iterator begin();
    std::array<double, 3>::iterator end();
    std::array<double, 3>::const_iterator cbegin() const;
    std::array<double, 3>::const_iterator cend() const;

    const double& operator[](size_t i) const noexcept;
    size_t size() const noexcept;
    const std::array<double, 3>& raw() const noexcept;

    Point& operator+=(const Point& other);
    Point& operator-=(const Point& other);
    Point& operator*=(double x);
    Point& operator/=(double x);
    Point operator+(const Point& other) const;
    Point operator-() const;
    Point operator-(const Point& other) const;
    Point operator*(double x) const;
    Point operator/(double x) const;

    double dist(const Point& other) const;
    bool isnear(const Point& other, double eps = 1e-5) const;
    bool isnan() const;

    friend bool operator==(const Point& lhs, const Point& rhs);
    friend bool operator!=(const Point& lhs, const Point& rhs);
    friend std::ostream& operator<<(std::ostream& os, const Point& p);

private:
    void copy();

    std::shared_ptr<VecImpl> pimpl_;
};

namespace std {
template <>
struct hash<Point> {
    size_t operator()(const Point& key) const
    {
        return hasher_(key.raw());
    }

private:
    hash<array<double, 3>> hasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_POINT_H_
