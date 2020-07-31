#include <gtest/gtest.h>

#include <gm/compare.hpp>
#include <gm/point.hpp>
#include <gm/surf_point.hpp>
#include <gm/surfaces.hpp>
#include <util/debug.hpp>
#include <util/math.hpp>

#include <random>

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

TEST(TestSurface, cylinder_projection_1)
{
    CylindricalSurface s(5, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {1, 1, 1}));

    static constexpr size_t niter = 1000;

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

TEST(TestSurface, toroidal_projection)
{
    ToroidalSurface s(0.5, 1, Axis::from_xy({1, 0, 0}, {0, 1, 0}, {0, 0, 0}));
    SurfPoint p0(1, 1);
    auto p = s.f(p0);
    auto q = s.project(p);
    ASSERT_NEAR(dist(p, s(q)), 0, 1e-5);
}
