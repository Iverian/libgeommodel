#include <geom_model/abstract_curve.h>
#include <geom_model/geom_util.h>
#include <util/debug.h>

#include <iostream>
#include <iterator>

#include <util/math.h>

#include <fmt/ostream.h>

using namespace std;

Vec AbstractCurve::df(double u) const
{
    return ::diff<Vec>([this](double t) { return Vec(f(t)); }, u);
}

Vec AbstractCurve::df2(double u) const
{
    return ::diff2<Vec>([this](double t) { return Vec(f(t)); }, u);
}

double AbstractCurve::length(double begin, double end) const
{
    static constexpr size_t max_iter = 10;
    static constexpr double eps = 1e-5;

    size_t i, mesh_size = 10;
    double prev, result = 0;
    for (i = 0; i < max_iter; ++i) {
        prev = result;
        result = approx_length(begin, end, mesh_size);
        mesh_size *= 2;
        if (fabs(prev - result) < eps) {
            break;
        }
    }
    if (i == max_iter) {
        fmt::print(cerr,
                   "[WARN]: ({0}): Iteration for curve {1} with parameters "
                   "{{ \"begin\": {2}, \"end\": {3} }} exceeded allowed limit",
                   __FUNCTION__, *this, begin, end);
    }
    return result;
}

optional<double> AbstractCurve::project_greater(const Point& x,
                                                double min) const
{
    optional<double> result = project(x);
    if (result.value() < min) {
        result = nullopt;
    }
    return result;
};

Point AbstractCurve::operator()(double u) const
{
    return f(u);
}

Vec AbstractCurve::tangent(double u) const
{
    return df(u);
}

Vec AbstractCurve::normal(double u) const
{
    auto d = df(u), d2 = df2(u);
    auto v = norm(d);
    auto w = ::diff<double>([this](double u) { return df(u).norm(); }, u);
    return (d2 - (w / v) * d) / sqr(v);
}

Point AbstractCurve::gproject(const Point& p) const
{
    return f(project(p));
}

double AbstractCurve::curvature(double u) const
{
    auto d = df(u);
    return cross(d, df2(u)).norm() / pow(d.norm(), 3);
}

optional<double> AbstractCurve::project_iterative(const Point& p, double init,
                                                  size_t max_iter) const
{
    static constexpr auto eps = 1e-5;

    optional<double> result = nullopt;
    double u = init, delta;
    for (size_t i = 0; i < max_iter; ++i) {
        auto v = p - f(u);
        auto d = df(u);
        delta = dot(d, v) / (dot(df2(u), v) - sqr(d));
        u -= delta;
        if (isnear(delta, 0, eps)) {
            result = u;
            break;
        }
    }
    return result;
}

bool AbstractCurve::is_init_in_interval(const Point& p,
                                        const array<double, 2>& interval) const
{
    array<bool, 2> flags;
    transform(cbegin(interval), cend(interval), begin(flags),
              [&](const auto& u) { return signbit(dot(p - f(u), df(u))); });
    return flags[0] != flags[1];
}

double AbstractCurve::approx_length(double begin, double end,
                                    size_t mesh_size) const
{
    auto step = (end - begin) / mesh_size;
    vector<double> y(mesh_size);
    for (size_t i = 0; i < mesh_size; ++i)
        y.emplace_back(norm(df(begin + i * step)));
    return ::trapz(y, step);
}

ostream& operator<<(ostream& os, const AbstractCurve& c)
{
    return c.print(os);
}
