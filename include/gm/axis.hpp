#ifndef GEOM_MODEL_INCLUDE_GM_AXIS_HPP_
#define GEOM_MODEL_INCLUDE_GM_AXIS_HPP_

#include "exports.h"
#include "point.h"
#include "vec.h"

#include <ostream>
#include <tuple>
#include <utility>

namespace gm {

class GM_EXPORT Axis {
    static constexpr size_t N = 3;

public:
    using View = std::tuple<Point, Vec, Vec, Vec>;

    Axis() noexcept;
    Axis(const Point& c, const std::array<Vec, N>& basis) noexcept;

    [[nodiscard]] const Vec& operator[](size_t i) const noexcept;
    [[nodiscard]] const Vec& x() const noexcept;
    [[nodiscard]] const Vec& y() const noexcept;
    [[nodiscard]] const Vec& z() const noexcept;
    [[nodiscard]] const Point& c() const noexcept;

    [[nodiscard]] View get_view() const noexcept;

    [[nodiscard]] Vec rotate_x(double angle, const Vec& v) const noexcept;
    [[nodiscard]] Vec rotate_y(double angle, const Vec& v) const noexcept;
    [[nodiscard]] Vec rotate_z(double angle, const Vec& v) const noexcept;

    [[nodiscard]] Point rotate_x(double angle, const Point& p) const noexcept;
    [[nodiscard]] Point rotate_y(double angle, const Point& p) const noexcept;
    [[nodiscard]] Point rotate_z(double angle, const Point& p) const noexcept;

    [[nodiscard]] Point global(const Point& p) const noexcept;
    [[nodiscard]] Vec global(const Vec& v) const noexcept;

    [[nodiscard]] Point pglobal(double x, double y, double z) const noexcept;
    [[nodiscard]] Vec vglobal(double x, double y, double z) const noexcept;

    friend bool operator==(const Axis& lhs, const Axis& rhs) noexcept;
    friend bool operator!=(const Axis& lhs, const Axis& rhs) noexcept;

    [[nodiscard]] static Axis from_xy(const Vec& x, const Vec& ref_y,
                                      const Point& c) noexcept;
    [[nodiscard]] static Axis from_zx(const Vec& z, const Vec& ref_x,
                                      const Point& c) noexcept;
    [[nodiscard]] static Axis from_abc(const Point& a, const Point& b,
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

#endif // GEOM_MODEL_INCLUDE_GM_AXIS_HPP_
