#include <geom_model/abstract_surface.h>
#include <geom_model/geom_util.h>
#include <geom_model/parametric_point.h>
#include <geom_model/plane.h>
#include <geom_model/point.h>
#include <geom_model/vec.h>

#include <util/math.h>

#include <algorithm>

using namespace std;

Point AbstractSurface::operator()(const ParametricPoint& p) const
{
    return f(p);
}

Vec AbstractSurface::normal(const ParametricPoint& p) const
{
    return cross(dfu(p), dfv(p));
}

Plane AbstractSurface::tangent(const ParametricPoint& p) const
{
    return Plane(Axis::from_xy(dfu(p), dfv(p), f(p)));
}

Point AbstractSurface::gproject(const Point& p) const
{
    return f(project(p));
}

function<Vec(double)> AbstractSurface::u_fixed(double v) const
{
    return [v, this](double u) { return Vec(f({u, v})); };
}

function<Vec(double)> AbstractSurface::v_fixed(double u) const
{
    return [u, this](double v) { return Vec(f({u, v})); };
}

Vec AbstractSurface::dfu(const ParametricPoint& p) const
{
    return ::diff<Vec>(u_fixed(p.v), p.u);
}

Vec AbstractSurface::dfv(const ParametricPoint& p) const
{
    return ::diff<Vec>(v_fixed(p.u), p.v);
}

Vec AbstractSurface::dfuu(const ParametricPoint& p) const
{
    return ::diff2<Vec>(u_fixed(p.v), p.u);
}

Vec AbstractSurface::dfvv(const ParametricPoint& p) const
{
    return ::diff2<Vec>(v_fixed(p.u), p.v);
}

Vec AbstractSurface::dfuv(const ParametricPoint& p) const
{
    return ::diff11<Vec>([this](ParametricPoint t) { return Vec(f(t)); }, p);
}

optional<ParametricPoint>
AbstractSurface::project_iterative(const Point& p, const ParametricPoint& init,
                                   size_t max_iter) const
{
    optional<ParametricPoint> result = nullopt;
    ParametricPoint u = init, d;
    for (size_t i = 0; i < max_iter; ++i) {
        d = project_to_step(p, u);
        u += d;
        if (isnear(max(abs(d)), 0, 1e-5)) {
            result = u;
            break;
        }
    }
    return result;
}

ParametricPoint
AbstractSurface::project_to_step(const Point& x,
                                 const ParametricPoint& p) const
{
    auto dv = Vec(x - f(p)), fu = dfu(p), fv = dfv(p);
    auto a11 = dot(dv, dfuu(p)) - sqr(fu),
         a12 = dot(dv, dfuv(p)) - dot(fu, fv),
         a22 = dot(dv, dfvv(p)) - sqr(fv);
    auto b1 = -dot(dv, fu), b2 = -dot(dv, fv);
    auto det = a11 * a22 - sqr(a12);
    auto det1 = b1 * a22 - b2 * a12;
    auto det2 = b2 * a11 - b1 * a12;

    return {det1 / det, det2 / det};
}

bool AbstractSurface::is_init_in_square(
    const Point& p, const array<ParametricPoint, 2>& diag) const
{
    auto c = (diag[0] + diag[1]) / 2;
    auto d = abs(diag[1] - diag[0]) / 2;

    auto fu = dfu(c), fv = dfv(c);
    auto a1 = dot(fu, p - f({c.u - d.u, c.v})),
         a2 = dot(fu, p - f({c.u + d.u, c.v})),
         b1 = dot(fv, p - f({c.u, c.v - d.v})),
         b2 = dot(fv, p - f({c.u, c.v + d.v}));
    auto a = a1 * a2, b = b1 * b2;
    return (a <= 0) && (b <= 0);
}

ostream& operator<<(ostream& os, const AbstractSurface& s)
{
    return s.print(os);
}
