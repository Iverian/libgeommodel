#include <gm/bspline_curve.h>
#include <gm/compare.h>
#include <gm/point.h>
#include <gm/vec.h>

#include <primitive/wpoint.h>
#include <util/math.h>
#include <util/util.h>
#include <util/debug.h>

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <vector>

#include <fmt/ostream.h>

#include "cox_de_boor.h"

using namespace std;

namespace gm {

struct BSplineCurve::Impl {
    Impl(size_t degree, vector<double> knots, vector<Point> points,
         vector<double> weights = vector<double>());
    Impl(size_t degree, const vector<size_t>& knot_mult,
         const vector<double>& knot_list, vector<Point> points,
         vector<double> weights = vector<double>());

    Point f(double u) const noexcept;
    Vec df(double u) const noexcept;
    Vec df2(double u) const noexcept;
    ostream& print(ostream& os) const;

    double project(const Point& p, const BSplineCurve& c) const;
    double ustep(double u) const;
    optional<double> newton_iter(const Point& p, double init,
                                 size_t max_iter) const;
    bool init_between(const Point& p, double a, double b) const;
    double bound_check(double u) const;

    double pfront() const noexcept;
    double pback() const noexcept;

private:
    void init_cpoints(const vector<Point>& p, const vector<double>& w);

    size_t order_;
    vector<double> knots_;
    vector<CPoint> cpoints_;
    CoxDeBoor<CPoint> cdb_;
};

BSplineCurve::~BSplineCurve() = default;

BSplineCurve::BSplineCurve(BSplineCurve&&) noexcept = default;

BSplineCurve& BSplineCurve::operator=(BSplineCurve&&) noexcept = default;

BSplineCurve::BSplineCurve(const BSplineCurve& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}

BSplineCurve& BSplineCurve::operator=(const BSplineCurve& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

BSplineCurve::BSplineCurve(size_t degree, vector<double> knots,
                           vector<Point> points, vector<double> weights)
    : pimpl_(make_unique<BSplineCurve::Impl>(degree, knots, points, weights))
{
}

BSplineCurve::BSplineCurve(size_t degree,
                           const vector<size_t>& knot_multiplies,
                           const vector<double>& knot_list,
                           vector<Point> points, vector<double> weights)
    : pimpl_(make_unique<BSplineCurve::Impl>(degree, knot_multiplies,
                                             knot_list, points))
{
}

Point BSplineCurve::f(double u) const noexcept
{
    return pimpl_->f(u);
}

Vec BSplineCurve::df(double u) const noexcept
{
    return pimpl_->df(u);
}

Vec BSplineCurve::df2(double u) const noexcept
{
    return pimpl_->df2(u);
}




ostream& BSplineCurve::print(ostream& os) const
{
    return pimpl_->print(os);
}

double BSplineCurve::project(const Point& p) const
{
    return pimpl_->project(p, *this);
}


double BSplineCurve::pfront() const noexcept {
    return pimpl_->pfront();
}

double BSplineCurve::pback() const noexcept {
    return pimpl_->pback();
}

BSplineCurve::Impl::Impl(size_t degree, vector<double> knots,
                         vector<Point> points, vector<double> weights)
    : order_(degree + 1)
    , knots_(move(knots))
    , cpoints_()
    , cdb_()
{
    init_cpoints(points, weights);
}

BSplineCurve::Impl::Impl(size_t degree, const vector<size_t>& knot_mult,
                         const vector<double>& knot_list, vector<Point> points,
                         vector<double> weights)
    : order_(degree + 1)
    , knots_()
    , cpoints_()
    , cdb_()
{
    for (size_t i = 0; i < knot_mult.size(); ++i) {
        for (size_t j = knot_mult[i]; j > 0; --j) {
            knots_.emplace_back(knot_list[i]);
        }
    }
    init_cpoints(points, weights);
}

Point BSplineCurve::Impl::f(double u) const noexcept
{
    return Point(cdb_.proxy(u).get(0).p());
}

Vec BSplineCurve::Impl::df(double u) const noexcept
{
    auto p = cdb_.proxy(u).range(2);
    return Vec(CPoint::d1(p).p()); // CPoint::d1(p).v();
}

Vec BSplineCurve::Impl::df2(double u) const noexcept
{
    auto p = cdb_.proxy(u).range(3);
    return Vec(CPoint::d2(p).p()); // CPoint::d2(p).v();
}

double BSplineCurve::Impl::pfront() const noexcept
{
    return knots_.front();
}

double BSplineCurve::Impl::pback() const noexcept
{
    return knots_.back();
}

ostream& BSplineCurve::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"deg\": {0}, \"knots\": {1}, "
               "\"cpoints\": {2} }}",
               order_ - 1, knots_, cpoints_);
    return os;
}

double BSplineCurve::Impl::project(const Point& p, const BSplineCurve& c) const
{
    static constexpr size_t dcoeff = 10;
    static constexpr size_t max_iter = 1000;

    optional<double> result;
    auto min_dist = numeric_limits<double>::max();

    auto ubegin = knots_.front();
    auto uend = knots_.back();
    auto udelta = (uend - ubegin) / dcoeff;

    auto u = ubegin;
    while (!isnear(u, uend)) {
        auto v = u + ::min(uend - u, udelta * ustep(u));
        if (init_between(p, u, v)) {
            for (auto& t : {u, v}) {
                auto op = newton_iter(p, t, max_iter);
                if (op.has_value()) {
                    auto q = op.value();
                    if (auto cur_dist = dist(f(q), p); cur_dist < min_dist) {
                        min_dist = cur_dist;
                        result = q;
                    }
                }
            }
        }
        u = v;
    }
    check_if(result.has_value(), "unable to project point {} on curve {}", p,
             c);

    return result.value();
}

bool BSplineCurve::Impl::init_between(const Point& p, double a, double b) const
{
    auto fa = signbit(dot(df(a), f(a) - p));
    auto fb = signbit(dot(df(b), f(b) - p));
    return fa != fb;
}

optional<double> BSplineCurve::Impl::newton_iter(const Point& p, double init,
                                                 size_t max_iter) const
{
    optional<double> result = nullopt;

    auto u = init;
    for (size_t i = 0; i < max_iter; ++i) {
        auto w = Vec(f(u) - p);
        auto d = df(u);
        auto v = bound_check(u - dot(d, w) / (dot(df2(u), w) + sqr(d)));

        if ((iszero(norm(w)) && iszero(cos(d, w)))
            || iszero(norm((v - u) * d))) {
            result = u;
            break;
        }

        u = v;
    }

    return result;
}

double BSplineCurve::Impl::bound_check(double u) const
{
    if (u > knots_.back()) {
        u = knots_.back();
    } else if (u < knots_.front()) {
        u = knots_.front();
    }
    return u;
}

double BSplineCurve::Impl::ustep(double u) const
{
    double cmul = 1;
    if (order_ > 2) {
        auto d1 = df(u);
        cmul = sqr(d1) / norm(cross(d1, df2(u)));
    }
    return cmul;
}

void BSplineCurve::Impl::init_cpoints(const vector<Point>& p,
                                      const vector<double>& w)
{
    auto n = p.size();
    cpoints_.resize(n);

    if (w.empty()) {
        for (size_t i = 0; i < n; ++i) {
            cpoints_[i] = CPoint(p[i].raw(), 1);
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            cpoints_[i] = CPoint(p[i].raw(), w[i]);
        }
    }

    cdb_ = CoxDeBoor<CPoint>(order_, knots_, cpoints_);
}

} // namespace gm