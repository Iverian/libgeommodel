#ifndef GEOM_MODEL_INCLUDE_GM_SURF_POINT_H_
#define GEOM_MODEL_INCLUDE_GM_SURF_POINT_H_

#include <algorithm>
#include <iostream>

namespace gm {

struct SurfPoint;

std::ostream& operator<<(std::ostream& os, const SurfPoint& p);
SurfPoint abs(const SurfPoint& p);

struct SurfPoint {
    constexpr SurfPoint(double pu, double pv)
        : u {pu}
        , v {pv}
    {
    }

    constexpr SurfPoint(const SurfPoint& other)
        : u {other.u}
        , v {other.v}
    {
    }

    constexpr SurfPoint(SurfPoint&& other) noexcept
        : u {other.u}
        , v {other.v}
    {
    }

    constexpr SurfPoint& operator=(const SurfPoint& other)
    {
        u = other.u;
        v = other.v;
        return *this;
    }

    constexpr SurfPoint& operator=(SurfPoint&& other) noexcept
    {
        u = std::move(other.u);
        v = std::move(other.v);
        return *this;
    }

    ~SurfPoint() = default;

    constexpr SurfPoint& operator+=(const SurfPoint& other)
    {
        u += other.u;
        v += other.v;
        return *this;
    }

    constexpr SurfPoint& operator-=(const SurfPoint& other)
    {
        u -= other.u;
        v -= other.v;
        return *this;
    }

    constexpr SurfPoint& operator*=(double x)
    {
        u *= x;
        v *= x;
        return *this;
    }

    constexpr SurfPoint& operator/=(double x)
    {
        u /= x;
        v /= x;
        return *this;
    }

    constexpr SurfPoint operator+(const SurfPoint& other) const
    {
        auto result = *this;
        return (result += other);
    }
    constexpr SurfPoint operator-(const SurfPoint& other) const
    {
        auto result = *this;
        return (result -= other);
    }

    constexpr SurfPoint operator*(double x) const
    {
        auto result = *this;
        return (result *= x);
    }

    constexpr SurfPoint operator/(double x) const
    {
        auto result = *this;
        return (result /= x);
    }

    constexpr SurfPoint()
        : u {0}
        , v {0}
    {
    }

    double u;
    double v;
};

constexpr SurfPoint operator*(double x, const SurfPoint& p)
{
    return (p * x);
}

constexpr bool operator==(const SurfPoint& lhs, const SurfPoint& rhs)
{
    return (lhs.u == lhs.u) && (rhs.v == rhs.v);
}

constexpr bool operator!=(const SurfPoint& lhs, const SurfPoint& rhs)
{
    return !(lhs == rhs);
}

constexpr double max(const SurfPoint& p)
{
    return std::max(p.u, p.v);
}

constexpr double min(const SurfPoint& p)
{
    return std::min(p.u, p.v);
}

} // namespace gm

#endif // GEOM_MODEL_INCLUDE_GM_SURF_POINT_H_