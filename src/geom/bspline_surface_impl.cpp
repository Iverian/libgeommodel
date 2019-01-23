#include "bspline_surface_impl.h"

#include <util/itertools.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

BSplineSurface::Impl::Impl()
    : order_()
    , knots_()
    , cpoints_()
    , cdb_()
{
}

BSplineSurface::Impl::Impl(const order_type& order, const knot_type& knots,
                           const cpoint_type& cpoints, const CPointSize& size)
    : order_(order)
    , knots_(knots)
    , cpoints_(cpoints)
    , size_(size)
    , cdb_()
{
    init_cdb();
}

BSplineSurface::Impl::Impl(const order_type& order, const knot_type& knots,
                           const std::vector<std::vector<CPoint>>& cpoints)
    : order_(order)
    , knots_(knots)
    , cpoints_()
    , cdb_()
{
    init_cpoints(cpoints);
}

BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<double>& ku,
                           const vector<double>& kv,
                           const std::vector<std::vector<Point>>& p,
                           const std::vector<std::vector<double>>& w)
{
    init_cpoints(w, p);
}

BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<size_t>& ku_mult,
                           const vector<double>& ku_vals,
                           const vector<size_t>& kv_mult,
                           const vector<double>& kv_vals,
                           const std::vector<std::vector<Point>>& p,
                           const std::vector<std::vector<double>>& w)
{
    for (size_t i = 0; i < ku_mult.size(); ++i)
        for (auto j = ku_mult[i]; j > 0; --j)
            knots_[0].emplace_back(ku_vals[i]);
    for (size_t i = 0; i < kv_mult.size(); ++i)
        for (auto j = kv_mult[i]; j > 0; --j)
            knots_[1].emplace_back(kv_vals[i]);
    for (auto& i : knots_)
        i.shrink_to_fit();

    init_cpoints(w, p);
}

Point BSplineSurface::Impl::f(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = cdb_[i].proxy(p.v).get(0);
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Point(cdb_u.proxy(p.u).get(0).p());
}

Vec BSplineSurface::Impl::dfu(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = cdb_[i].proxy(p.v).get(0);
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Vec(CPoint::d1(cdb_u.proxy(p.u).range(2)).p());
}

Vec BSplineSurface::Impl::dfv(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = CPoint::d1(cdb_[i].proxy(p.v).range(2));
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Vec(cdb_u.proxy(p.u).get(0).p());
}

Vec BSplineSurface::Impl::dfuu(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = cdb_[i].proxy(p.v).get(0);
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Vec(CPoint::d2(cdb_u.proxy(p.u).range(3)).p());
}

Vec BSplineSurface::Impl::dfvv(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = CPoint::d2(cdb_[i].proxy(p.v).range(3));
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Vec(cdb_u.proxy(p.u).get(0).p());
}

Vec BSplineSurface::Impl::dfuv(const SurfPoint& p) const noexcept
{
    auto n = size_.n;
    vector<CPoint> cp(n);

    for (size_t i = 0; i < n; ++i) {
        cp[i] = CPoint::d1(cdb_[i].proxy(p.v).range(2));
    }

    auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
    return Vec(CPoint::d1(cdb_u.proxy(p.u).range(2)).p());
}

ostream& BSplineSurface::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"du\": {}, \"dv\": {}, \"ku\": "
               "{}, \"kv\": {}, \"cpoints\": {}, \"shape\": [{}, {}] }}",
               order_[0] - 1, order_[1] - 1,
               RangePrint(begin(knots_[0]), end(knots_[0])),
               RangePrint(begin(knots_[1]), end(knots_[1])),
               RangePrint(begin(cpoints_), end(cpoints_)), size_.n, size_.m);
    return os;
}

// TODO: реализовать
SurfPoint BSplineSurface::Impl::project(const Point& p) const
{
    return SurfPoint();
}

const BSplineSurface::Impl::order_type& BSplineSurface::Impl::order() const
    noexcept
{
    return order_;
}

const BSplineSurface::Impl::knot_type& BSplineSurface::Impl::knots() const
    noexcept
{
    return knots_;
}

const BSplineSurface::Impl::cpoint_type& BSplineSurface::Impl::cpoints() const
    noexcept
{
    return cpoints_;
}

const BSplineSurface::Impl::CPointSize& BSplineSurface::Impl::size() const
    noexcept
{
    return size_;
}

BSplineSurface::Impl BSplineSurface::Impl::to_bezier() const
{
    Impl result;
    auto count = order_;
    auto prev = knots_[0];

    result.order_ = order_;
    return *this;
}

#define _(i, j) ((size_.n) * (i) + (j))

void BSplineSurface::Impl::init_cpoints(
    const std::vector<std::vector<CPoint>>& cp)
{
    size_.n = cp.size();
    size_.m = cp[0].size();
    cpoints_.resize(size_.n * size_.m);

    for (size_t i = 0; i < size_.n; ++i) {
        for (size_t j = 0; j < size_.m; ++j) {
            cpoints_[_(i, j)] = cp[i][j];
        }
    }

    init_cdb();
}

void BSplineSurface::Impl::init_cpoints(
    const std::vector<std::vector<double>>& w,
    const std::vector<std::vector<Point>>& p)
{
    size_.n = p.size();
    size_.m = p[0].size();
    cpoints_.resize(size_.n * size_.m);

    if (w.empty()) {
        for (size_t i = 0; i < size_.n; ++i) {
            for (size_t j = 0; j < size_.m; ++j) {
                cpoints_[_(i, j)] = CPoint(p[i][j].raw(), 1);
            }
        }
    } else {
        for (size_t i = 0; i < size_.n; ++i) {
            for (size_t j = 0; j < size_.m; ++j) {
                cpoints_[_(i, j)] = CPoint(p[i][j].raw(), w[i][j]);
            }
        }
    }

    init_cdb();
}

#undef _

void BSplineSurface::Impl::init_cdb()
{
    cdb_.resize(size_.n);

    VectorView<double> knots(knots_[1]);
    auto data = cpoints_.data();
    for (size_t i = 0; i < size_.n; ++i) {
        cdb_[i] = CoxDeBoor<CPoint>(order_[1], knots,
                                    VectorView<CPoint>(data, size_.m));
        data += size_.m;
    }
}

} // namespace gm