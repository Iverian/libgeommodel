#ifndef GEOM_MODEL_INCLUDE_AXIS_H_
#define GEOM_MODEL_INCLUDE_AXIS_H_

#include "mat.h"

#include <array>
#include <iostream>
#include <memory>
#include <tuple>

class Axis;

std::ostream& operator<<(std::ostream& os, const Axis& ax);

class Axis {
public:
    ~Axis();
    Axis(Axis&&) noexcept;
    Axis& operator=(Axis&&) noexcept;
    Axis(const Axis& other);
    Axis& operator=(const Axis& other);

    Axis();
    Axis(const Vec& x, const Vec& y, const Vec& z, const Point& center);

    const Vec& at(size_t j) const;
    const Vec& operator[](size_t j) const;
    const Vec& x() const;
    const Vec& y() const;
    const Vec& z() const;
    const Point& center() const;
    const Mat& get_transform() const;
    std::tuple<Point, Vec, Vec, Vec> get_view() const;
    double arg(const Point& p) const;

    Vec rotate_x(double angle, const Vec& v);
    Vec rotate_y(double angle, const Vec& v);
    Vec rotate_z(double angle, const Vec& v);

    Point rotate_x(double angle, const Point& p);
    Point rotate_y(double angle, const Point& p);
    Point rotate_z(double angle, const Point& p);

    Point global(const Point& p) const;
    Vec global(const Vec& v) const;
    Point pglobal(double x, double y, double z) const;
    Vec vglobal(double x, double y, double z) const;

    friend bool operator==(const Axis& lhs, const Axis& rhs);
    friend bool operator!=(const Axis& lhs, const Axis& rhs);

    static Axis from_xy(const Vec& x, const Vec& ref_y, const Point& center);
    static Axis from_zx(const Vec& z, const Vec& ref_x, const Point& center);
    static Axis from_abc(const Point& a, const Point& b, const Point& c);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

namespace std {
template <>
struct hash<Axis> {
    size_t operator()(const Axis& key) const
    {
        return (mhasher_(key.get_transform()) << 1) ^ phasher_(key.center());
    }

private:
    hash<Mat> mhasher_;
    hash<Point> phasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_AXIS_H_
