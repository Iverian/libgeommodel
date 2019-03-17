#include <bspline/bspline_curve_impl.hpp>
#include <gm/bspline_curve.hpp>

using namespace std;

namespace gm {

BSplineCurve::BSplineCurve(size_t degree, ::vector<double> knots,
                           ::vector<Point> points, ::vector<double> weights)
    : pimpl_(make_unique<BSplineCurve::Impl>(degree, knots, points, weights))
{
}

BSplineCurve::BSplineCurve(size_t degree,
                           const ::vector<size_t>& knot_multiplies,
                           const ::vector<double>& knot_list,
                           ::vector<Point> points, ::vector<double> weights)
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
    return pimpl_->project(p);
}

double BSplineCurve::pfront() const noexcept
{
    return pimpl_->pfront();
}

double BSplineCurve::pback() const noexcept
{
    return pimpl_->pback();
}

} // namespace gm