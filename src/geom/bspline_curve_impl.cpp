#include "bspline_curve_impl.h"
#include "curve_projector.h"

#include <util/itertools.h>

#include <fmt/ostream.h>

using namespace std;

namespace gm {

BSplineCurve::Impl::~Impl() = default;

BSplineCurve::Impl::Impl(BSplineCurve::Impl&&) noexcept = default;

BSplineCurve::Impl& BSplineCurve::Impl::
operator=(BSplineCurve::Impl&&) noexcept
    = default;

BSplineCurve::Impl::Impl(const BSplineCurve::Impl& rhs)
    : order_(rhs.order_)
    , knots_(rhs.knots_)
    , cpoints_(rhs.cpoints_)
    , cdb_()
    , proj_(rhs.proj_ ? make_unique<CurveProjector>(*rhs.proj_) : nullptr)
{
    init_cdb();
}

BSplineCurve::Impl& BSplineCurve::Impl::
operator=(const BSplineCurve::Impl& rhs)
{
    order_ = rhs.order_;
    knots_ = rhs.knots_;
    cpoints_ = rhs.cpoints_;
    proj_ = rhs.proj_ ? make_unique<CurveProjector>(*rhs.proj_) : nullptr;

    init_cdb();

    return *this;
};

BSplineCurve::Impl::Impl(const BezierPatch& patch)
    : order_(patch.cp.size())
    , knots_(2 * order_)
    , cpoints_(patch.cp)
    , cdb_()
    , proj_(nullptr)
{
    for (size_t i = 0; i < order_; ++i) {
        knots_[i] = patch.front;
        knots_[order_ + i] = patch.back;
    }

    init_cdb();
}

BSplineCurve::Impl::Impl(size_t degree, vector<double> knots,
                         vector<Point> points, vector<double> weights)
    : order_(degree + 1)
    , knots_(move(knots))
    , cpoints_()
    , cdb_()
    , proj_(nullptr)
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
    , proj_(nullptr)
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
    return Vec(CPoint::d1(p).p());
}

Vec BSplineCurve::Impl::df2(double u) const noexcept
{
    auto p = cdb_.proxy(u).range(3);
    return Vec(CPoint::d2(p).p());
}

double BSplineCurve::Impl::pfront() const noexcept
{
    return knots_.front();
}

double BSplineCurve::Impl::pback() const noexcept
{
    return knots_.back();
}

// см. NURBS book 5.3 knot refinement
vector<BezierPatch> BSplineCurve::Impl::bezier_patches() const
{
    vector<BezierPatch> result;
    vector<double> alphas(order_);

    auto p = order_ - 1;        // степень кривой
    auto m = knots_.size() - 1; // cpoints_.size() + p;

    auto a = p;
    auto b = p + 1;
    // bool completed = false;

    size_t nb = 0;

    result.emplace_back(order_);
    for (size_t i = 0; i < order_; ++i) {
        result[0].cp[i] = cpoints_[i];
    }

    while (b < m) {
        auto i = b;
        while (b < m && isnear(knots_[b], knots_[b + 1], Tolerance::MAX)) {
            ++b;
        }
        result[nb].front = knots_[a];
        result[nb].back = knots_[b];

        if (b < m) {
            result.emplace_back(order_);
        } else {
            // completed = true;
        }

        auto mult = b - i + 1;
        if (mult < p) {
            auto numer = knots_[b] - knots_[a];
            for (size_t j = p; j > mult; --j) {
                alphas[j - mult - 1] = numer / (knots_[a + j] - knots_[a]);
            }

            auto r = p - mult;
            for (size_t j = 1; j < r + 1; ++j) {
                auto save = r - j;
                auto s = mult + j;
                for (size_t k = p; k >= s; --k) {
                    auto& alpha = alphas[k - s];
                    result[nb].cp[k] = alpha * result[nb].cp[k]
                        + (1. - alpha) * result[nb].cp[k - 1];
                }

                if (b < m) {
                    result[nb + 1].cp[save] = result[nb].cp[p];
                }
            }
            ::fill(::begin(alphas), ::end(alphas), 0.);
        }

        ++nb;
        if (b < m) {
            for (size_t i = p - mult; i < order_; ++i) {
                result[nb].cp[i] = cpoints_[b - p + i];
            }
            a = b;
            ++b;
        }
    }

    return result;
}

ostream& BSplineCurve::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"deg\": {0}, \"knots\": {1}, "
               "\"cpoints\": {2} }}",
               order_ - 1, RangePrint(begin(knots_), end(knots_)),
               RangePrint(begin(cpoints_), end(cpoints_)));
    return os;
}

double BSplineCurve::Impl::project(const Point& p) const
{
    if (!proj_) {
        proj_.reset(new CurveProjector(*this));
    }

    return proj_->get(p);
    // static constexpr size_t dcoeff = 10;
    // static constexpr size_t max_iter = 1000;

    // optional<double> result;
    // auto min_dist = numeric_limits<double>::max();

    // auto ubegin = knots_.front();
    // auto uend = knots_.back();
    // auto udelta = (uend - ubegin) / dcoeff;

    // auto u = ubegin;
    // while (!isnear(u, uend)) {
    //     auto v = u + ::min(uend - u, udelta * ustep(u));
    //     if (init_between(p, u, v)) {
    //         for (auto& t : {u, v}) {
    //             auto op = newton_iter(p, t, max_iter);
    //             if (op.has_value()) {
    //                 auto q = op.value();
    //                 if (auto cur_dist = dist(f(q), p); cur_dist < min_dist)
    //                 {
    //                     min_dist = cur_dist;
    //                     result = q;
    //                 }
    //             }
    //         }
    //     }
    //     u = v;
    // }
    // check_if(result.has_value(), "unable to project point {} on curve {}",
    // p,
    //          c);

    // return result.value();
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

    init_cdb();
}

void BSplineCurve::Impl::init_cdb()
{
    cdb_ = CoxDeBoor<CPoint>(order_, knots_, cpoints_);
}

BezierPatch::BezierPatch(size_t order)
    : front(0)
    , back(1)
    , cp(order)
{
}

size_t BezierPatch::order() const noexcept
{
    return cp.size();
}

} // namespace gm