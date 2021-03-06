#include <gm/abstract_curve.hpp>
#include <gm/compare.hpp>

#include <iostream>
#include <iterator>

#include <util/math.hpp>

#include <fmt/ostream.h>

namespace gm {

AbstractCurve::~AbstractCurve() = default;

Vec AbstractCurve::df(double u) const noexcept
{
    return ::diff<Vec>([this](double t) { return Vec(f(t)); }, u);
}

Vec AbstractCurve::df2(double u) const noexcept
{
    return ::diff2<Vec>([this](double t) { return Vec(f(t)); }, u);
}

std::optional<double> AbstractCurve::project_greater(const Point& x,
                                                     double min) const noexcept
{
    std::optional<double> result = project(x);
    if (result.value() < min) {
        result = std::nullopt;
    }
    return result;
};

Point AbstractCurve::operator()(double u) const noexcept
{
    return f(u);
}

Vec AbstractCurve::tangent(double u) const noexcept
{
    return df(u);
}

Vec AbstractCurve::normal(double u) const noexcept
{
    auto d = df(u), d2 = df2(u);
    auto v = norm(d);
    auto w = diff<double>([this](double u) { return norm(df(u)); }, u);
    return (d2 - (w / v) * d) / sqr(v);
}

Vec AbstractCurve::unit_normal(double u) const __GM_NOEXCEPT_RELEASE__
{
    return unit(normal(u));
}

Point AbstractCurve::gproject(const Point& p) const noexcept
{
    return f(project(p));
}

double AbstractCurve::curvature(double u) const noexcept
{
    auto d = df(u);
    return norm(cross(d, df2(u))) / pow(norm(d), 3);
}

double AbstractCurve::approx_length(double begin, double end, size_t n) const
{
    auto step = (end - begin) / n;
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i)
        y[i] = norm(df(begin + i * step));
    return trapz(std::begin(y), std::end(y), step);
}

std::ostream& operator<<(std::ostream& os, const AbstractCurve& c)
{
    return c.print(os);
}

} // namespace gm