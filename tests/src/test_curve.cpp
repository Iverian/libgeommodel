#include <gtest/gtest.h>

#include <gm/curves.hpp>
#include <gm/point.hpp>
#include <util/debug.hpp>

#include <random>
#include <stdexcept>

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
    Line c(Vec(10000, 100, 0), Point(1, 2, 3));

    std::default_random_engine rnd;
    std::uniform_real_distribution<double> udist;
    Vec z(0, 0, 1);

    for (size_t i = 0; i < 100; ++i) {
        auto t = udist(rnd);

        auto p = c.f(t);
        auto x = p + (2. + udist(rnd)) * z;
        auto r = c.project(x);
        auto q = c.f(r);
        auto d = c.df(r);

        debug_fmt(std::cout, "TEST 1: POINT #{} t: {}", i, t);

        auto qx = Vec(q, x);
        ASSERT_LE(dist(q, x), dist(p, x));
        ASSERT_NEAR(cos(qx, d), 0., cmp::tol());
    }
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
    std::cout << "t=" << unit(t) << ", n=" << unit(n) << std::endl;

    ASSERT_NEAR(dot(t, n), 0, 1e-5);
}
