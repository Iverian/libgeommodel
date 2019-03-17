#define _USE_MATH_DEFINES

#include <gtest/gtest.h>

#include <gm/curves.hpp>
#include <gm/surfaces.hpp>
#include <util/debug.hpp>

#include <cmath>

using namespace std;
using namespace gm;

TEST(TestPrimitive, vec_operations)
{
    Vec u(1, 2, 3), v(1, 0, 3);

    ASSERT_EQ(u + v, Vec(2, 2, 6));
    ASSERT_EQ(u - v, Vec(0, 2, 0));
    ASSERT_EQ(2 * u, Vec(2, 4, 6));
    ASSERT_EQ((-3) * v, Vec(-3, 0, -9));
    ASSERT_EQ(0 * u, Vec(0, 0, 0));
    ASSERT_EQ(u / 3, Vec(1. / 3, 2. / 3, 1));
    ASSERT_EQ(dot(u, v), 10);
    ASSERT_EQ(cross(u, v), Vec(6, 0, -2));
    ASSERT_EQ(norm(u), sqrt(14));
    ASSERT_EQ(norm(v), sqrt(10));
}

TEST(TestPrimitive, point_operations)
{
    Point u(1, 2, 3), v(1, 0, 3);

    ASSERT_EQ(u + v, Point(2, 2, 6));
    ASSERT_EQ(u - v, Point(0, 2, 0));
    ASSERT_EQ(2 * u, Point(2, 4, 6));
    ASSERT_EQ((-3) * v, Point(-3, 0, -9));
    ASSERT_EQ(0 * u, Point(0, 0, 0));
    ASSERT_EQ(u / 3, Point(1. / 3, 2. / 3, 1));
    ASSERT_EQ(dot(u, v), 10);
    ASSERT_EQ(dist(u, v), 2);
}

TEST(TestPrimitive, diff_correct)
{
    Parabola c;
    auto d = [](double u) { return 2 * (u * Vec(1, 0, 0) + Vec(0, 1, 0)); };
    ASSERT_NEAR(norm(c.df(3) - d(3)), 0, 1e-5);
    ASSERT_NEAR(norm(c.df2(3) - Vec(2, 0, 0)), 0, 1e-5);
}

TEST(TestPrimitive, diff2_correct)
{
    double r = 1000;
    Vec x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);
    Point c(0, 0, 0);
    SphericalSurface s(r, Axis::from_xy(x, y, c));
    SurfPoint p(2, 1);

    ASSERT_NEAR(norm(s.dfu(p) - s.AbstractSurface::dfu(p)), 0, 1e-2);
    ASSERT_NEAR(norm(s.dfv(p) - s.AbstractSurface::dfv(p)), 0, 1e-2);
    ASSERT_NEAR(norm(s.dfuu(p) - s.AbstractSurface::dfuu(p)), 0, 1e-2);
    ASSERT_NEAR(norm(s.dfuv(p) - s.AbstractSurface::dfuv(p)), 0, 1e-2);
    ASSERT_NEAR(norm(s.dfvv(p) - s.AbstractSurface::dfvv(p)), 0, 1e-2);
}

TEST(TestPrimitive, init_axis)
{
    Vec z(0, -1, 0);
    Vec x(1, 0, 0);
    auto y = cross(z, x);
    Point c(1, 1, 1);
    auto ax = Axis::from_zx(z, x, c);

    ASSERT_EQ(ax[0], x);
    ASSERT_EQ(ax[1], y);
    ASSERT_EQ(ax[2], z);
    ASSERT_EQ(ax.c(), c);
}

TEST(TestPrimitive, rotation)
{
    auto ax = Axis::from_zx({0, 0, 1}, {1, 0, 0}, {0, 0, 0});
    auto a = unit({1, 0, 0});
    auto b = unit({1, 1, 0});
    auto c = M_PI_2;
    debug_fmt(cout, "r(a): {}, r(b): {}", ax.rotate_z(c, a),
              ax.rotate_z(c, b));

    SUCCEED();
}