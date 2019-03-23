#include <gtest/gtest.h>

#include <gm/compare.hpp>
#include <gm/point.hpp>

#include <util/debug.hpp>
#include <util/math.hpp>

#include <bspline/bspline_curve_impl.hpp>
#include <bspline/bspline_surface_impl.hpp>
#include <bspline/curve_projector.hpp>
#include <bspline/distance_surface.hpp>
#include <gm/surf_point.hpp>

#include <random>

using namespace gm;

static constexpr size_t N = 100;

struct TestBSplineImpl : ::testing::Test {
protected:
    TestBSplineImpl()
        : c(3, {4, 1, 1, 1, 1, 1, 1, 1, 4},
            {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.},
            {Point(100, 10, 0), Point(90, 10, -50), Point(80, 10, 0),
             Point(70, 10, 30), Point(60, 10, 0), Point(50, 10, -30),
             Point(40, 10, 0), Point(30, 10, 20), Point(20, 10, 0),
             Point(10, 10, -20), Point(0, 10, 0)})
        , s()
        , bz(c.bezier_patches())
        , bzs()

    {
        std::array<Point, 4> cp = {Point(0.99131233, -0.04928702, 0),
                                   Point(1.07925827, 1.32291813, 0),
                                   Point(-1.07925827, 1.32291813, 0),
                                   Point(-0.99131233, -0.04928702, 0)};
        Vec z(0, 0, 1);
        s = BSplineSurface::Impl(3, 1, {4, 4}, {0, 1}, {2, 2}, {0, 1},
                                 {{cp[0] + z, cp[0]},
                                  {cp[1] + z, cp[1]},
                                  {cp[2] + z, cp[2]},
                                  {cp[3] + z, cp[3]}});
        bzs = s.bezier_patches();
    }

    BSplineCurve::Impl c;
    BSplineSurface::Impl s;
    std::vector<BSplineCurve::Impl::BezierPatch> bz;
    std::vector<BSplineSurface::Impl::BezierPatch> bzs;
};

TEST_F(TestBSplineImpl, test_bspline_patches)
{
    for (auto& patch : bz) {
        auto impl = BSplineCurve::Impl(patch);
        auto step = (patch.pback() - patch.pfront()) / (N - 1);
        auto u = patch.pfront();

        for (size_t i = 0; i < N; ++i) {
            EXPECT_EQ(impl.f(u), c.f(u));
            EXPECT_EQ(impl.df(u), c.df(u));
            // EXPECT_EQ(impl.df2(u), c.df2(u));

            u += step;
        }
    }

    SUCCEED();
}

TEST_F(TestBSplineImpl, test_distance_curve)
{
    Point p(1, 1, 1);
    auto f = [this, &p](double u) { return sqr(c.f(u) - p); };

    for (auto& patch : bz) {
        auto g = DistanceCurve(patch, p);
        auto step = (patch.pback() - patch.pfront()) / (N - 1);
        auto u = patch.pfront();

        for (size_t i = 0; i < N; ++i) {
            EXPECT_NEAR(f(u), g.f(u), cmp::tol());

            u += step;
        }
    }
}

TEST_F(TestBSplineImpl, test_distance_surface)
{
    Point p(100, 100, 100);
    auto f = [this, &p](SurfPoint r) { return sqr(s.f(r) - p); };

    for (auto& patch : bzs) {
        auto g = DistanceSurface(patch, p);
        auto h = (patch.pback() - patch.pfront()) / (N - 1);
        auto& q = patch.pfront();

        auto r = q;
        for (size_t i = 0; i < N; ++i) {
            r.v = q.v;
            for (size_t j = 0; j < N; ++j) {
                EXPECT_NEAR(f(r), g.f(r), cmp::tol());
                r.v += h.v;
            }
            r.u += h.u;
        }
    }
}
