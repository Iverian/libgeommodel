#include <gtest/gtest.h>

#include <gm/bspline_curve.h>
#include <gm/bspline_surface.h>
#include <gm/compare.h>
#include <gm/point.h>
#include <util/debug.h>
#include <util/math.h>

#include <random>

using namespace std;
using namespace gm;

double fpad(double t, double eps) noexcept;

class TestBSpline : public ::testing::Test {
protected:
    static constexpr size_t niter = 1000;
    static constexpr auto pad = 1e-2;

    TestBSpline()
        : c(3, {4, 1, 1, 1, 1, 1, 1, 1, 4},
            {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.},
            {Point(100, 10, 0), Point(90, 10, -50), Point(80, 10, 0),
             Point(70, 10, 30), Point(60, 10, 0), Point(50, 10, -30),
             Point(40, 10, 0), Point(30, 10, 20), Point(20, 10, 0),
             Point(10, 10, -20), Point(0, 10, 0)})
        , s()
    {
        array<Point, 4> cp = {Point(0.99131233, -0.04928702, 0),
                              Point(1.07925827, 1.32291813, 0),
                              Point(-1.07925827, 1.32291813, 0),
                              Point(-0.99131233, -0.04928702, 0)};
        Vec z(0, 0, 1);
        s = BSplineSurface(3, 1, {4, 4}, {0, 1}, {2, 2}, {0, 1},
                           {{cp[0] + z, cp[0]},
                            {cp[1] + z, cp[1]},
                            {cp[2] + z, cp[2]},
                            {cp[3] + z, cp[3]}});
    }

    BSplineCurve c;
    BSplineSurface s;
};

TEST_F(TestBSpline, eval)
{
    static constexpr auto eps = 1e-1;
    static constexpr auto eps2 = 1e-3;

    auto& f = s;

    ASSERT_NEAR(dist(f({0, 0}), Point(1, 0, 1)), 0, eps);
    ASSERT_NEAR(dist(f({0, 1}), Point(1, 0, 0)), 0, eps);
    ASSERT_NEAR(dist(f({0.5, 1}), Point(0, 1, 0)), 0, eps);

    array<SurfPoint, 6> test_points
        = {SurfPoint(0.1, 0.1), SurfPoint(0.1, 0.9), SurfPoint(0.9, 0.1),
           SurfPoint(0.5, 0.5), SurfPoint(0.1, 0.5), SurfPoint(0.5, 0.1)};

    for (auto& p : test_points) {
        auto a = f.dfu(p);
        auto b = f.AbstractSurface::dfu(p);
        // debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfv(p);
        auto b = f.AbstractSurface::dfv(p);
        // debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfuu(p);
        auto b = f.AbstractSurface::dfuu(p);
        // debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfuv(p);
        auto b = f.AbstractSurface::dfuv(p);
        // debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfvv(p);
        auto b = f.AbstractSurface::dfvv(p);
        // debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }

    SUCCEED();
}

TEST_F(TestBSpline, curve_proj_1)
{
    default_random_engine rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = udist(rnd);

        debug_fmt(cout, "TEST 1: POINT #{} t: {}", i, t);
        auto p = c.f(t);
        auto r = c.project(p);
        auto q = c.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, curve_proj_2)
{
    default_random_engine rnd;
    uniform_real_distribution<double> udist;
    Vec z(0, 0, 1);

    for (size_t i = 0; i < niter; ++i) {
        auto t = fpad(udist(rnd), pad);

        debug_fmt(cout, "TEST 2: POINT #{} t: {}", i, t);
        auto p = c.f(t);
        auto x = p + (2. + udist(rnd)) * z;
        auto r = c.project(x);
        auto q = c.f(r);
        auto d = c.df(r);

        auto qx = Vec(q, x);
        debug_fmt(cout, "r: {} qx: {}", r, qx);
        ASSERT_NEAR(cos(qx, d), 0., cmp::tol());
    }
}

TEST_F(TestBSpline, surface_proj_1)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(fpad(udist(rnd), pad), fpad(udist(rnd), pad));

        debug_fmt(cout, "TEST 1: POINT #{} t: {}", i, t);
        auto p = s.f(t);
        auto r = s.project(p);
        auto q = s.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

#undef PAD_

TEST_F(TestBSpline, surface_proj_2)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(rnd() % 2, udist(rnd));

        debug_fmt(cout, "TEST 2: POINT #{} t: {}", i, t);
        auto p = s.f(t);
        auto r = s.project(p);
        auto q = s.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, surface_proj_3)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(udist(rnd), udist(rnd));
        auto p = s.f(t);
        auto x = p + (0.5 + udist(rnd)) * s.unit_normal(t);

        debug_fmt(cout, "TEST 3: POINT #{} t: {} p: {} p + a*n: {}", i, t, p,
                  x);
        auto r = s.project(x);
        auto q = s.f(r);

        auto d_xq = dist(x, q);
        auto d_xp = dist(x, p);
        debug_fmt(cout, "r: {} dist(x, q): {} dist(x, p): {}", r, d_xq, d_xp);
        ASSERT_TRUE(d_xq < d_xp || cmp::near(d_xp, d_xq, Tolerance::SINGLE));
    }
}

inline double fpad(double t, double eps) noexcept
{
    return (1 - 2 * eps) * t + eps;
}