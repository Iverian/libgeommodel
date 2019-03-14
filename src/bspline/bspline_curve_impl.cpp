#include <bspline/bspline_curve_impl.h>
#include <bspline/curve_projector.h>
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
    : c_(rhs.c_)
    , proj_(rhs.proj_ ? make_unique<CurveProjector>(*rhs.proj_) : nullptr)
{
}

BSplineCurve::Impl& BSplineCurve::Impl::
operator=(const BSplineCurve::Impl& rhs)
{
    c_ = rhs.c_;
    proj_ = rhs.proj_ ? make_unique<CurveProjector>(*rhs.proj_) : nullptr;

    return *this;
};

BSplineCurve::Impl::Impl(const BezierPatch& patch)
    : c_(patch.order(), patch.knots(), patch.cpoints())
    , proj_(nullptr)
{
}

BSplineCurve::Impl::Impl(size_t degree, vector<double> knots,
                         vector<Point> points, vector<double> weights)
    : c_()
    , proj_(nullptr)
{
    init_curve(degree + 1, knots, points, weights);
}

BSplineCurve::Impl::Impl(size_t degree, const vector<size_t>& knot_mult,
                         const vector<double>& knot_list, vector<Point> points,
                         vector<double> weights)
    : c_()
    , proj_(nullptr)
{
    vector<double> k;
    for (size_t i = 0; i < knot_mult.size(); ++i) {
        for (size_t j = knot_mult[i]; j > 0; --j) {
            k.emplace_back(knot_list[i]);
        }
    }

    init_curve(degree + 1, k, points, weights);
}

Point BSplineCurve::Impl::f(double u) const noexcept
{
    return Point(c_.f(u));
}

Vec BSplineCurve::Impl::df(double u) const noexcept
{
    return Vec(c_.df(u));
}

Vec BSplineCurve::Impl::df2(double u) const noexcept
{
    return Vec(c_.df2(u));
}

double BSplineCurve::Impl::pfront() const noexcept
{
    return c_.pfront();
}

double BSplineCurve::Impl::pback() const noexcept
{
    return c_.pback();
}

vector<BSplineCurve::Impl::BezierPatch>
BSplineCurve::Impl::bezier_patches() const
{
    return c_.bezier_patches();
}

ostream& BSplineCurve::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"deg\": {0}, \"knots\": {1}, "
               "\"cpoints\": {2} }}",
               c_.order() - 1, RangePrint(begin(c_.knots()), end(c_.knots())),
               RangePrint(begin(c_.cpoints()), end(c_.cpoints())));
    return os;
}

double BSplineCurve::Impl::project(const Point& p) const
{
    if (!proj_) {
        proj_.reset(new CurveProjector(*this));
    }

    return proj_->call(p);
}

void BSplineCurve::Impl::init_curve(size_t order, const vector<double>& k,
                                    const vector<Point>& p,
                                    const vector<double>& w)
{
    auto n = p.size();
    vector<Super::CPoint> cp(n);

    if (w.empty()) {
        for (size_t i = 0; i < n; ++i) {
            cp[i] = Super::CPoint(p[i].raw(), 1);
        }
    } else {
        for (size_t i = 0; i < n; ++i) {
            cp[i] = Super::CPoint(p[i].raw(), w[i]);
        }
    }

    c_ = Super(order, k, cp);
}

} // namespace gm