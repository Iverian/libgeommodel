#include <gtest/gtest.h>

#include <gm/bspline_curve.hpp>
#include <gm/bspline_surface.hpp>
#include <gm/compare.hpp>
#include <gm/point.hpp>
#include <util/debug.hpp>
#include <util/math.hpp>

#include <random>

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
        , c1(3, {0, 0, 0, 0, 0.5, 1, 1, 1, 1},
             {gm::Point(-10, 40, 10), gm::Point(-10, 10, 10),
              gm::Point(0, -10, 10), gm::Point(10, 10, 10),
              gm::Point(10, 40, 10)})
        , s()
    {
        std::array<Point, 4> cp = {Point(0.99131233, -0.04928702, 0),
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
    BSplineCurve c1;
    BSplineSurface s;
};

TEST_F(TestBSpline, print)
{
    EXPECT_NO_THROW({ fmt::print(std::cout, "{}", s); });
}

TEST_F(TestBSpline, curve_diff)
{
    static constexpr auto eps = 1e-2;
    auto& f = c1;

    for (double u = 0.; gm::cmp::le(u, 1.); u += eps) {
        auto b = f.df2(u);
        debug_fmt(std::cout, "{} : {}", u, b);
    }
}

TEST_F(TestBSpline, eval_surface)
{
    static constexpr auto eps = 1e-1;
    static constexpr auto eps2 = 1e-3;

    auto& f = s;

    ASSERT_NEAR(dist(f({0, 0}), Point(1, 0, 1)), 0, eps);
    ASSERT_NEAR(dist(f({0, 1}), Point(1, 0, 0)), 0, eps);
    ASSERT_NEAR(dist(f({0.5, 1}), Point(0, 1, 0)), 0, eps);

    fmt::print(std::cout, "{}", s);
    std::array<SurfPoint, 6> test_points
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
    std::default_random_engine rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = udist(rnd);

        debug_fmt(std::cout, "TEST 1: POINT #{} t: {}", i, t);
        auto p = c.f(t);
        auto r = c.project(p);
        auto q = c.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(std::cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, curve_proj_2)
{
    std::default_random_engine rnd;
    std::uniform_real_distribution<double> udist;
    Vec z(0, 0, 1);

    for (size_t i = 0; i < niter; ++i) {
        auto t = fpad(udist(rnd), pad);

        debug_fmt(std::cout, "TEST 2: POINT #{} t: {}", i, t);
        auto p = c.f(t);
        auto x = p + (2. + udist(rnd)) * z;
        auto r = c.project(x);
        auto q = c.f(r);
        auto d = c.df(r);

        auto qx = Vec(q, x);
        debug_fmt(std::cout, "r: {} qx: {}", r, qx);
        ASSERT_LE(dist(q, x), dist(p, x));
        ASSERT_NEAR(cos(qx, d), 0., cmp::tol());
    }
}

TEST_F(TestBSpline, surface_proj_1)
{
    std::minstd_rand0 rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(fpad(udist(rnd), pad), fpad(udist(rnd), pad));

        debug_fmt(std::cout, "TEST 1: POINT #{} t: {}", i, t);
        auto p = s.f(t);
        auto r = s.project(p);
        auto q = s.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(std::cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, surface_proj_2)
{
    std::minstd_rand0 rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(rnd() % 2, udist(rnd));

        debug_fmt(std::cout, "TEST 2: POINT #{} t: {}", i, t);
        auto p = s.f(t);
        auto r = s.project(p);
        auto q = s.f(r);

        auto d_pq = dist(p, q);
        debug_fmt(std::cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(cmp::zero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, surface_proj_3)
{
    std::minstd_rand0 rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(udist(rnd), udist(rnd));
        auto p = s.f(t);
        auto x = p + udist(rnd) * s.unit_normal(t);

        debug_fmt(std::cout, "TEST 3: POINT #{} t: {} p: {} p + a*n: {}", i, t,
                  p, x);
        EXPECT_NO_THROW({
            auto r = s.project(x);
            auto q = s.f(r);
            auto aqx = angle(Vec(q, x), s.normal(r));
            auto dqx = dist(q, x);
            auto dpx = dist(p, x);

            debug_fmt(std::cout,
                      " r: {} dist(p, q): {} dist(q, x): {} angle(qx): {}", r,
                      dist(p, q), dist(q, x), aqx);
            EXPECT_PRED3(cmp::le, dqx, dpx, cmp::default_tolerance);
            EXPECT_NEAR(aqx, 0., cmp::tol(Tolerance::SINGLE));
        });
    }
}

TEST_F(TestBSpline, surface_proj_4)
{
    std::minstd_rand0 rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(rnd() % 2, udist(rnd));
        auto p = s.f(t);
        auto x = p + udist(rnd) * s.unit_normal(t);

        debug_fmt(std::cout,
                  "TEST 3: POINT #{} t: {} p: {} p + a*n: {} dist(p, x): {}",
                  i, t, p, x, dist(p, x));
        EXPECT_NO_THROW({
            auto r = s.project(x);
            auto q = s.f(r);
            auto aqx = angle(Vec(q, x), s.normal(r));

            debug_fmt(std::cout,
                      " r: {} dist(p, q): {} dist(q, x): {} angle(qx): {}", r,
                      dist(p, q), dist(q, x), aqx);
            EXPECT_TRUE(cmp::le(dist(q, x), dist(p, x)));
            EXPECT_NEAR(aqx, 0., cmp::tol(Tolerance::SINGLE));
        });
    }
}

TEST_F(TestBSpline, surface_proj_5)
{
    std::minstd_rand0 rnd;
    std::uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(udist(rnd), rnd() % 2);
        auto p = s.f(t);
        auto x = p + udist(rnd) * s.unit_normal(t);

        debug_fmt(std::cout,
                  "TEST 3: POINT #{} t: {} p: {} p + a*n: {} dist(p, x): {}",
                  i, t, p, x, dist(p, x));
        EXPECT_NO_THROW({
            auto r = s.project(x);
            auto q = s.f(r);
            auto aqx = angle(Vec(q, x), s.normal(r));

            debug_fmt(std::cout,
                      " r: {} dist(p, q): {} dist(q, x): {} angle(qx): {}", r,
                      dist(p, q), dist(q, x), aqx);
            EXPECT_TRUE(cmp::le(dist(q, x), dist(p, x)));
            EXPECT_NEAR(aqx, 0., cmp::tol(Tolerance::SINGLE));
        });
    }
}

inline double fpad(double t, double eps) noexcept
{
    return (1 - 2 * eps) * t + eps;
}