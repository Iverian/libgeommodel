#include <geom_model/bspline_surface.h>

#include "cox_de_boor.h"

#include <util/debug.h>
#include <util/itertools.h>
#include <util/math.h>
#include <util/util.h>

#include <algorithm>
#include <stdexcept>

#include <fmt/ostream.h>

using namespace std;

struct BSplineSurface::Impl {
    Impl(size_t du, size_t dv, const vector<double>& ku,
         const vector<double>& kv, const vector<vector<Point>>& p,
         const vector<vector<double>>& w = vector<vector<double>>());

    Impl(size_t du, size_t dv, const vector<size_t>& ku_mult,
         const vector<double>& ku_vals, const vector<size_t>& kv_mult,
         const vector<double>& kv_vals, const vector<vector<Point>>& p,
         const vector<vector<double>>& w = vector<vector<double>>());

    Point f(const ParametricPoint& p) const;
    Vec dfu(const ParametricPoint& p) const;
    Vec dfv(const ParametricPoint& p) const;
    Vec dfuu(const ParametricPoint& p) const;
    Vec dfvv(const ParametricPoint& p) const;
    Vec dfuv(const ParametricPoint& p) const;

    ostream& print(ostream& os) const;
    ParametricPoint project(const Point& p,
                            const BSplineSurface& parent) const;
    void init_weights();
    void init_cdb();

private:
    array<size_t, 2> order_;
    array<vector<double>, 2> knots_;
    vector<vector<Point>> control_points_;
    vector<vector<double>> weights_;
    vector<CoxDeBoor> cdb_;
};

BSplineSurface::BSplineSurface(size_t du, size_t dv, const vector<double>& ku,
                               const vector<double>& kv,
                               const vector<vector<Point>>& p,
                               const vector<vector<double>>& w)
    : pimpl_(make_unique<BSplineSurface::Impl>(du, dv, ku, kv, p, w))
{
}

BSplineSurface::BSplineSurface(size_t du, size_t dv,
                               const vector<size_t>& ku_mult,
                               const vector<double>& ku_vals,
                               const vector<size_t>& kv_mult,
                               const vector<double>& kv_vals,
                               const vector<vector<Point>>& p,
                               const vector<vector<double>>& w)
    : pimpl_(make_unique<BSplineSurface::Impl>(du, dv, ku_mult, ku_vals,
                                               kv_mult, kv_vals, p, w))
{
}

Point BSplineSurface::f(const ParametricPoint& p) const
{
    return pimpl_->f(p);
}

Vec BSplineSurface::dfu(const ParametricPoint& p) const
{
    return pimpl_->dfu(p);
}

Vec BSplineSurface::dfv(const ParametricPoint& p) const
{
    return pimpl_->dfv(p);
}

Vec BSplineSurface::dfuu(const ParametricPoint& p) const
{
    return pimpl_->dfuu(p);
}

Vec BSplineSurface::dfuv(const ParametricPoint& p) const
{
    return pimpl_->dfuv(p);
}

Vec BSplineSurface::dfvv(const ParametricPoint& p) const
{
    return pimpl_->dfvv(p);
}

ostream& BSplineSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

ParametricPoint BSplineSurface::project(const Point& p) const
{
    return pimpl_->project(p, *this);
}

BSplineSurface::~BSplineSurface() = default;
BSplineSurface::BSplineSurface(BSplineSurface&&) noexcept = default;
BSplineSurface& BSplineSurface::operator=(BSplineSurface&&) noexcept = default;
BSplineSurface::BSplineSurface(const BSplineSurface& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_)){};
BSplineSurface& BSplineSurface::operator=(const BSplineSurface& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
};

BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<double>& ku,
                           const vector<double>& kv,
                           const vector<vector<Point>>& p,
                           const vector<vector<double>>& w)
    : order_{du + 1, dv + 1}
    , knots_{ku, kv}
    , control_points_(p)
    , weights_(w)
    , cdb_()
{
    init_weights();
    init_cdb();
}

BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<size_t>& ku_mult,
                           const vector<double>& ku_vals,
                           const vector<size_t>& kv_mult,
                           const vector<double>& kv_vals,
                           const vector<vector<Point>>& p,
                           const vector<vector<double>>& w)
    : order_{du + 1, dv + 1}
    , knots_()
    , control_points_(p)
    , weights_(w)
    , cdb_()
{
    for (size_t i = 0; i < ku_mult.size(); ++i)
        for (auto j = ku_mult[i]; j > 0; --j)
            knots_[0].emplace_back(ku_vals[i]);
    for (size_t i = 0; i < kv_mult.size(); ++i)
        for (auto j = kv_mult[i]; j > 0; --j)
            knots_[1].emplace_back(kv_vals[i]);
    for (auto& i : knots_)
        i.shrink_to_fit();

    init_weights();
    init_cdb();
}

Point BSplineSurface::Impl::f(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n);

    for (size_t i = 0; i < n; ++i) {
        auto q = cdb_[i].get_proxy(p.v).get(0);
        r[i] = q.r / q.w;
        w[i] = q.w;
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return cdb_u.get_proxy(p.u).f();
}

Vec BSplineSurface::Impl::dfu(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n);

    for (size_t i = 0; i < n; ++i) {
        auto q = cdb_[i].get_proxy(p.v).get(0);
        r[i] = q.r / q.w;
        w[i] = q.w;
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return cdb_u.get_proxy(p.u).df();
}

Vec BSplineSurface::Impl::dfv(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n, 1);

    for (size_t i = 0; i < n; ++i) {
        r[i] = Point(cdb_[i].get_proxy(p.v).df());
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return Vec(cdb_u.get_proxy(p.u).f());
}

Vec BSplineSurface::Impl::dfuu(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n);

    for (size_t i = 0; i < n; ++i) {
        auto q = cdb_[i].get_proxy(p.v).get(0);
        r[i] = q.r / q.w;
        w[i] = q.w;
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return cdb_u.get_proxy(p.u).df2();
}

Vec BSplineSurface::Impl::dfvv(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n, 1);

    for (size_t i = 0; i < n; ++i) {
        r[i] = Point(cdb_[i].get_proxy(p.v).df2());
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return Vec(cdb_u.get_proxy(p.u).f());
}

Vec BSplineSurface::Impl::dfuv(const ParametricPoint& p) const
{
    auto n = control_points_.size();
    vector<Point> r(n);
    vector<double> w(n, 1);

    for (size_t i = 0; i < n; ++i) {
        r[i] = Point(cdb_[i].get_proxy(p.v).df());
    }

    auto cdb_u = CoxDeBoor(order_[0], knots_[0], r, w);
    return cdb_u.get_proxy(p.u).df();
}

ostream& BSplineSurface::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"du\": {0}, \"dv\": {1}, \"ku\": "
               "{2}, \"kv\": {3}, \"points\": {4}, \"weights\": {5} }}",
               order_[0] - 1, order_[1] - 1, knots_[0], knots_[1],
               control_points_, weights_);
    return os;
}

#define mesh_element(type, i)                                                 \
    (knots_[type].front() + (knots_[type].back() * (i)) / mesh_size[type])

ParametricPoint
BSplineSurface::Impl::project(const Point& p,
                              const BSplineSurface& parent) const
{
    static constexpr array<size_t, 2> mesh_size = {10, 10};

    vector<ParametricPoint> results;
    for (size_t i = 0; i < mesh_size[0]; ++i) {
        auto a = ParametricPoint{mesh_element(0, i), 0},
             b = ParametricPoint{mesh_element(0, i + 1), 0};
        for (size_t j = 0; j < mesh_size[1]; ++j) {
            a.v = mesh_element(1, j);
            b.v = mesh_element(1, j + 1);
            if (parent.is_init_in_square(p, {a, b})) {
                auto projection = parent.project_iterative(p, (a + b) / 2);
                if (projection.has_value()) {
                    results.emplace_back(projection.value());
                }
            }
        }
    }
    if (results.empty()) {
        throw runtime_error(fmt::format(
            "Unable to project point {0} to surface {1}", p, parent));
    }
    return *min_element(cbegin(results), cend(results),
                        [&](const auto& lhs, const auto& rhs) {
                            return dist(f(lhs), p) < dist(f(rhs), p);
                        });
}

#undef mesh_element

void BSplineSurface::Impl::init_weights()
{
    if (weights_.empty()) {
        auto n = control_points_.size();
        auto m = control_points_.front().size();
        weights_.resize(n);
        for (size_t i = 0; i < n; ++i)
            weights_[i].resize(m, 1);
    }
}

void BSplineSurface::Impl::init_cdb()
{
    cdb_.resize(control_points_.size());
    for (size_t i = 0; i < cdb_.size(); ++i) {
        cdb_[i]
            = CoxDeBoor(order_[1], knots_[1], control_points_[i], weights_[i]);
    }
}