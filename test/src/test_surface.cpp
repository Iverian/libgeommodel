#include <gtest/gtest.h>

#include <geom_model/point.h>
#include <geom_model/surfaces.h>
#include <util/debug.h>
#include <util/math.h>

using namespace std;

TEST(TestSurface, b_spline_init)
{
    static constexpr auto eps = 1e-1;
    static constexpr auto eps2 = 1e-3;

    array<Point, 4> cp = {
        Point(0.99131233, -0.04928702, 0), Point(1.07925827, 1.32291813, 0),
        Point(-1.07925827, 1.32291813, 0), Point(-0.99131233, -0.04928702, 0)};
    Vec z(0, 0, 1);
    auto f = BSplineSurface(3, 1, {4, 4}, {0, 1}, {2, 2}, {0, 1},
                            {{cp[0] + z, cp[0]},
                             {cp[1] + z, cp[1]},
                             {cp[2] + z, cp[2]},
                             {cp[3] + z, cp[3]}});

    ASSERT_NEAR(dist(f({0, 0}), Point(1, 0, 1)), 0, eps);
    ASSERT_NEAR(dist(f({0, 1}), Point(1, 0, 0)), 0, eps);
    ASSERT_NEAR(dist(f({0.5, 1}), Point(0, 1, 0)), 0, eps);

    array<ParametricPoint, 6> test_points
        = {ParametricPoint(0.1, 0.1), ParametricPoint(0.1, 0.9),
           ParametricPoint(0.9, 0.1), ParametricPoint(0.5, 0.5),
           ParametricPoint(0.1, 0.5), ParametricPoint(0.5, 0.1)};

    for (auto& p : test_points) {
        ASSERT_NEAR(dist(f.dfu(p), f.AbstractSurface::dfu(p)), 0, eps2);
    }
    for (auto& p : test_points) {
        ASSERT_NEAR(dist(f.dfv(p), f.AbstractSurface::dfv(p)), 0, eps2);
    }
    for (auto& p : test_points) {
        ASSERT_NEAR(dist(f.dfuu(p), f.AbstractSurface::dfuu(p)), 0, eps2);
    }
    for (auto& p : test_points) {
        ASSERT_NEAR(dist(f.dfuv(p), f.AbstractSurface::dfuv(p)), 0, eps2);
    }
    for (auto& p : test_points) {
        ASSERT_NEAR(dist(f.dfvv(p), f.AbstractSurface::dfvv(p)), 0, eps2);
    }
}

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
    ParametricPoint p0(0, 1);
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
    ParametricPoint p0(1, 1);
    auto p = s.f(p0);
    auto q = s.project(p);
    ASSERT_NEAR(dist(p, s(q)), 0, 1e-5);
}

TEST(TestSurface, bspline_projection)
{
    array<Point, 4> cp = {
        Point(0.99131233, -0.04928702, 0), Point(1.07925827, 1.32291813, 0),
        Point(-1.07925827, 1.32291813, 0), Point(-0.99131233, -0.04928702, 0)};
    Vec z(0, 0, 1);
    auto s = BSplineSurface(3, 1, {4, 4}, {0, 1}, {2, 2}, {0, 1},
                            {{cp[0] + z, cp[0]},
                             {cp[1] + z, cp[1]},
                             {cp[2] + z, cp[2]},
                             {cp[3] + z, cp[3]}});
    ParametricPoint p0(0.5, 0.5);
    auto p = s.f(p0);
    auto q = s.project(p);
    ASSERT_NEAR(dist(p, s(q)), 0, 1e-5);
}