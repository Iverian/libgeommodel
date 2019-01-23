#ifndef GEOM_MODEL_INCLUDE_AXIS_H_
#define GEOM_MODEL_INCLUDE_AXIS_H_

#include "point.h"
#include "vec.h"

#include <ostream>
#include <tuple>
#include <utility>

namespace gm {

class Axis {
    static constexpr size_t N = 3;

public:
    using View = std::tuple<Point, Vec, Vec, Vec>;

    Axis() noexcept;
    Axis(const Point& c, const std::array<Vec, N>& basis) noexcept;

    const Vec& operator[](size_t i) const noexcept;
    const Vec& x() const noexcept;
    const Vec& y() const noexcept;
    const Vec& z() const noexcept;
    const Point& c() const noexcept;

    View get_view() const noexcept;

    Vec rotate_x(double angle, const Vec& v) const noexcept;
    Vec rotate_y(double angle, const Vec& v) const noexcept;
    Vec rotate_z(double angle, const Vec& v) const noexcept;

    Point rotate_x(double angle, const Point& p) const noexcept;
    Point rotate_y(double angle, const Point& p) const noexcept;
    Point rotate_z(double angle, const Point& p) const noexcept;

    Point global(const Point& p) const noexcept;
    Vec global(const Vec& v) const noexcept;

    Point pglobal(double x, double y, double z) const noexcept;
    Vec vglobal(double x, double y, double z) const noexcept;

    friend bool operator==(const Axis& lhs, const Axis& rhs) noexcept;
    friend bool operator!=(const Axis& lhs, const Axis& rhs) noexcept;

    static Axis from_xy(const Vec& x, const Vec& ref_y,
                        const Point& c) noexcept;
    static Axis from_zx(const Vec& z, const Vec& ref_x,
                        const Point& c) noexcept;
    static Axis from_abc(const Point& a, const Point& b,
                         const Point& c) noexcept;

private:
    Point c_;
    Vec basis_[N];
};

std::ostream& operator<<(std::ostream& os, const Axis& ax);

} // namespace gm

namespace std {
template <>
struct hash<gm::Axis> {
    size_t operator()(const gm::Axis& key) const
    {
        return (vhasher_(key.x()) << 3) ^ (vhasher_(key.y()) << 2)
            ^ (vhasher_(key.z()) << 1) ^ phasher_(key.c());
    }

private:
    hash<gm::Vec> vhasher_;
    hash<gm::Point> phasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_AXIS_H_
