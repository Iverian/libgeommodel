#include <gtest/gtest.hpp>

#include <gm/curves.hpp>

#include <stdexcept>

using namespace std;
using namespace gm;

TEST(TestCurve, circle_project)
{
    Point p(0, 2, 0);
    Ellipse curve;
    auto u = curve.project(p);
    ASSERT_DOUBLE_EQ(dist(p, curve(u)), 1);
}

TEST(TestCurve, line_project)
{
    auto p = Point(2, 1, 0);
    Line curve;
    auto u = curve.project(p);
    ASSERT_DOUBLE_EQ(dist(p, curve(u)), 1);
}

TEST(TestCurve, b_spline_init)
{
    auto f = BSplineCurve(3, {4, 4}, {0, 1},
                          {Point(0.99131233, -0.04928702, 0),
                           Point(1.07925827, 1.32291813, 0),
                           Point(-1.07925827, 1.32291813, 0),
                           Point(-0.99131233, -0.04928702, 0)});

    ASSERT_NEAR(dist(f(0), Point(0.99131233, -0.04928702, 0)), 0, 1e-3);
    ASSERT_NEAR(dist(f(0.5), Point(4.44089210e-16, 9.79866841e-01, 0)), 0,
                1e-3);
    ASSERT_NEAR(dist(f(1), Point(-0.99131233, -0.04928702, 0)), 0, 1e-3);

    ASSERT_NEAR(dist(f.df(0), Vec(0.26383782, 4.11661544, 0)), 0, 1e-3);
    ASSERT_NEAR(dist(f.df(0.5), Vec(-3.10585589e+00, -3.55271368e-15, 0)), 0,
                1e-3);
    ASSERT_NEAR(dist(f.df(1), Vec(0.26383782, -4.11661544, 0)), 0, 1e-3);

    ASSERT_NEAR(dist(f.df2(0), Vec(-13.47877485, -8.23323088, 0)), 0, 1e-3);
    ASSERT_NEAR(dist(f.df2(0.5), Vec(8.88178420e-16, -8.23323088e+00, 0)), 0,
                1e-3);
    ASSERT_NEAR(dist(f.df2(1), Vec(13.47877485, -8.23323088, 0)), 0, 1e-3);
}

TEST(TestCurve, normal_test)
{
    auto ax = Axis::from_zx(Vec(0, 0, 1), Vec(1, 0, 0), Point(0, 0, 0));
    Ellipse curve(1, 2, ax);
    double u = 1;
    auto t = curve.tangent(u), n = curve.normal(u);
    cout << "t=" << unit(t) << ", n=" << unit(n) << endl;

    ASSERT_NEAR(dot(t, n), 0, 1e-5);
}
