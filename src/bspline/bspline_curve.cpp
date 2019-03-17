#include <bspline/bspline_curve_impl.hpp>
#include <gm/bspline_curve.hpp>

namespace gm {

BSplineCurve::BSplineCurve(size_t degree, std::vector<double> knots,
                           std::vector<Point> points,
                           std::vector<double> weights)
    : pimpl_(
        std::make_unique<BSplineCurve::Impl>(degree, knots, points, weights))
{
}

BSplineCurve::BSplineCurve(size_t degree,
                           const std::vector<size_t>& knot_multiplies,
                           const std::vector<double>& knot_list,
                           std::vector<Point> points,
                           std::vector<double> weights)
    : pimpl_(std::make_unique<BSplineCurve::Impl>(degree, knot_multiplies,
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

std::ostream& BSplineCurve::print(std::ostream& os) const
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