#include <gtest/gtest.h>

#include <gm/bspline_surface.h>
#include <gm/compare.h>
#include <gm/point.h>
#include <util/debug.h>
#include <util/math.h>

#include <random>

using namespace std;
using namespace gm;

class TestBSpline : public ::testing::Test {
protected:
    static constexpr size_t niter = 1000;
    static constexpr auto pad = 1e-9;

    void SetUp()
    {
        array<Point, 4> cp = {Point(0.99131233, -0.04928702, 0),
                              Point(1.07925827, 1.32291813, 0),
                              Point(-1.07925827, 1.32291813, 0),
                              Point(-0.99131233, -0.04928702, 0)};
        Vec z(0, 0, 1);
        auto c = BSplineSurface(3, 1, {4, 4}, {0, 1}, {2, 2}, {0, 1},
                                {{cp[0] + z, cp[0]},
                                 {cp[1] + z, cp[1]},
                                 {cp[2] + z, cp[2]},
                                 {cp[3] + z, cp[3]}});
        c_ = make_shared<BSplineSurface>(move(c));
    }

    void TearDown()
    {
        c_.reset();
    }

    shared_ptr<BSplineSurface> c_;
};

TEST_F(TestBSpline, eval)
{
    static constexpr auto eps = 1e-1;
    static constexpr auto eps2 = 1e-3;

    auto& f = *c_;

    ASSERT_NEAR(dist(f({0, 0}), Point(1, 0, 1)), 0, eps);
    ASSERT_NEAR(dist(f({0, 1}), Point(1, 0, 0)), 0, eps);
    ASSERT_NEAR(dist(f({0.5, 1}), Point(0, 1, 0)), 0, eps);

    array<SurfPoint, 6> test_points
        = {SurfPoint(0.1, 0.1), SurfPoint(0.1, 0.9), SurfPoint(0.9, 0.1),
           SurfPoint(0.5, 0.5), SurfPoint(0.1, 0.5), SurfPoint(0.5, 0.1)};

    for (auto& p : test_points) {
        auto a = f.dfu(p);
        auto b = f.AbstractSurface::dfu(p);
        debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfv(p);
        auto b = f.AbstractSurface::dfv(p);
        debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfuu(p);
        auto b = f.AbstractSurface::dfuu(p);
        debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfuv(p);
        auto b = f.AbstractSurface::dfuv(p);
        debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
    for (auto& p : test_points) {
        auto a = f.dfvv(p);
        auto b = f.AbstractSurface::dfvv(p);
        debug_fmt(cout, "a: {} b: {}", a, b);
        ASSERT_NEAR(dist(a, b), 0, eps2);
    }
}

#define PAD_(x) ((1 + 2 * pad) * (x) + pad)

TEST_F(TestBSpline, projection_1)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(PAD_(udist(rnd)), PAD_(udist(rnd)));

        debug_fmt(cout, "TEST 1: POINT #{} t: {}", i, t);
        auto p = c_->f(t);
        auto r = c_->project(p);
        auto q = c_->f(r);

        auto d_pq = dist(p, q);
        debug_fmt(cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(iszero(d_pq, Tolerance::SINGLE));
    }
}

#undef PAD_

TEST_F(TestBSpline, projection_2)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(rnd() % 2, udist(rnd));

        debug_fmt(cout, "TEST 2: POINT #{} t: {}", i, t);
        auto p = c_->f(t);
        auto r = c_->project(p);
        auto q = c_->f(r);

        auto d_pq = dist(p, q);
        debug_fmt(cout, "r: {} dist(p, q): {}", r, d_pq);
        ASSERT_TRUE(iszero(d_pq, Tolerance::SINGLE));
    }
}

TEST_F(TestBSpline, projection_3)
{
    minstd_rand0 rnd;
    uniform_real_distribution<double> udist;

    for (size_t i = 0; i < niter; ++i) {
        auto t = SurfPoint(udist(rnd), udist(rnd));
        auto p = c_->f(t);
        auto x = p + (0.5 + udist(rnd)) * c_->unit_normal(t);

        debug_fmt(cout, "TEST 3: POINT #{} t: {} p: {} p + a*n: {}", i, t, p, x);
        auto r = c_->project(x);
        auto q = c_->f(r);

        auto d_xq = dist(x, q);
        auto d_xp = dist(x, p);
        debug_fmt(cout, "r: {} dist(x, q): {} dist(x, p): {}", r, d_xq, d_xp);
        ASSERT_TRUE(d_xq < d_xp || isnear(d_xp, d_xq, Tolerance::SINGLE));
    }
}
