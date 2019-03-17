#include <gtest/gtest.h>

#include <gm/compare.hpp>
#include <gm/point.hpp>
#include <gm/surfaces.hpp>
#include <util/debug.hpp>
#include <util/math.hpp>

#include <random>

using namespace std;
using namespace gm;

TEST(TestSurface, plane_projection)
{
    Plane s(Axis::from_xy({0, 1, 0}, {0, 0, 1}, {0, 0, 0}));
    Point p(5, 4, 2);
    auto q = s.project(p);
    ASSERT_NEAR(dist(s(q), Point(0, 4, 2)), 0, 1e-5);
}

TEST(TestSurface, sphere_projection)
{
    SphericalSurface s(1, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
    SurfPoint p0(0, 1);
    Point pp(1, 1, 1);
    auto p = s.f(p0);
    auto q = s.project(p);
    ASSERT_NEAR(dist(s(q), p), 0, 1e-5);
}

TEST(TestSurface, sphere_projection_1)
{
    SphericalSurface s(1, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
    Point p(1, 1, 1);
    auto q = s.gproject(p);
    double v = 1. / sqrt(3);
    ASSERT_NEAR(dist(q, Point(v, v, v)), 0, 1e-5);
}

TEST(TestSurface, conic_projection)
{
    ConicalSurface s(1, M_PI / 6,
                     Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
}

TEST(TestSurface, cylinder_projection)
{
    CylindricalSurface s(1, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
    Point p(1, 1, 1);
    auto q = s.project(p);
    auto v = 1 / sqrt(2);
    ASSERT_NEAR(dist(s(q), Point(v, v, 1)), 0, 1e-5);
}

TEST(TestSurface, toroidal_projection)
{
    ToroidalSurface s(0.5, 1, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
    SurfPoint p0(1, 1);
    auto p = s.f(p0);
    auto q = s.project(p);
    ASSERT_NEAR(dist(p, s(q)), 0, 1e-5);
}
