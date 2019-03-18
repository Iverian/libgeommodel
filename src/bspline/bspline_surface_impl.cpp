#include <bspline/bspline_surface_impl.hpp>
#include <bspline/surface_projector.hpp>
#include <bspline/wpoint.hpp>
#include <gm/bspline_surface.hpp>
#include <util/itertools.hpp>

#include <fmt/ostream.h>

namespace gm {

BSplineSurface::Impl::Impl()
    : c_()
    , proj_(nullptr)
{
}

BSplineSurface::Impl::~Impl() = default;
BSplineSurface::Impl::Impl(Impl&&) noexcept = default;
BSplineSurface::Impl& BSplineSurface::Impl::operator=(Impl&&) noexcept
    = default;
BSplineSurface::Impl::Impl(const Impl& rhs)
    : c_(rhs.c_)
    , proj_(rhs.proj_ ? std::make_unique<SurfaceProjector>(*rhs.proj_)
                      : nullptr)
{
}
BSplineSurface::Impl& BSplineSurface::Impl::operator=(const Impl& rhs)
{
    c_ = rhs.c_;
    proj_
        = rhs.proj_ ? std::make_unique<SurfaceProjector>(*rhs.proj_) : nullptr;

    return *this;
}

BSplineSurface::Impl::Impl(size_t du, size_t dv, const std::vector<double>& ku,
                           const std::vector<double>& kv,
                           const std::vector<std::vector<Point>>& p,
                           const std::vector<std::vector<double>>& w)
    : c_()
    , proj_(nullptr)
{
    KnotsType k {ku, kv};
    init_surface({du + 1, dv + 1}, std::move(k), p, w);
}

BSplineSurface::Impl::Impl(size_t du, size_t dv,
                           const std::vector<size_t>& ku_mult,
                           const std::vector<double>& ku_vals,
                           const std::vector<size_t>& kv_mult,
                           const std::vector<double>& kv_vals,
                           const std::vector<std::vector<Point>>& p,
                           const std::vector<std::vector<double>>& w)
    : c_()
    , proj_(nullptr)
{
    KnotsType k;

    for (size_t i = 0; i < ku_mult.size(); ++i) {
        for (auto j = ku_mult[i]; j > 0; --j) {
            k.first.emplace_back(ku_vals[i]);
        }
    }
    k.first.shrink_to_fit();
    for (size_t i = 0; i < kv_mult.size(); ++i) {
        for (auto j = kv_mult[i]; j > 0; --j) {
            k.second.emplace_back(kv_vals[i]);
        }
    }
    k.second.shrink_to_fit();

    init_surface({du + 1, dv + 1}, std::move(k), p, w);
}

Point BSplineSurface::Impl::f(const SurfPoint& p) const noexcept
{
    return Point(c_.f(p));
}

Vec BSplineSurface::Impl::dfu(const SurfPoint& p) const noexcept
{
    return Vec(c_.dfu(p));
}

Vec BSplineSurface::Impl::dfv(const SurfPoint& p) const noexcept
{
    return Vec(c_.dfv(p));
}

Vec BSplineSurface::Impl::dfuu(const SurfPoint& p) const noexcept
{
    return Vec(c_.dfuu(p));
}

Vec BSplineSurface::Impl::dfvv(const SurfPoint& p) const noexcept
{
    return Vec(c_.dfvv(p));
}

Vec BSplineSurface::Impl::dfuv(const SurfPoint& p) const noexcept
{
    return Vec(c_.dfuv(p));
}

std::ostream& BSplineSurface::Impl::print(std::ostream& os) const
{
    auto s = c_.shape();
    fmt::print(
        os,
        "{{ \"type\": \"bspline\", \"du\": {}, \"dv\": {}, \"ku\": "
        "{}, \"kv\": {}, \"cpoints\": {}, \"shape\": [{}, {}] }}",
        c_.order().first - 1, c_.order().second - 1,
        RangePrint(std::begin(c_.knots().first), std::end(c_.knots().first)),
        RangePrint(std::begin(c_.knots().first), std::end(c_.knots().second)),
        RangePrint(std::begin(c_.cpoints()), std::end(c_.cpoints())), s.first,
        s.second);
    return os;
}

SurfPoint BSplineSurface::Impl::project(const Point& p) const
{
    if (!proj_) {
        proj_.reset(new SurfaceProjector(*this));
    }

    return proj_->call(p);
}

const BSplineSurface::Impl::OrderType& BSplineSurface::Impl::order() const
    noexcept
{
    return c_.order();
}

const BSplineSurface::Impl::KnotsType& BSplineSurface::Impl::knots() const
    noexcept
{
    return c_.knots();
}

const BSplineSurface::Impl::CPointsType& BSplineSurface::Impl::cpoints() const
    noexcept
{
    return c_.cpoints();
}

std::vector<BSplineSurface::Impl::BezierPatch>
BSplineSurface::Impl::bezier_patches() const noexcept
{
    return c_.bezier_patches();
}

#define _w(i, j) ((w.empty()) ? (1.) : (w[(i)][(j)]))

void BSplineSurface::Impl::init_surface(
    std::pair<size_t, size_t> order,
    const std::pair<std::vector<double>, std::vector<double>>& k,
    const std::vector<std::vector<Point>>& p,
    const std::vector<std::vector<double>>& w)
{
    std::pair<size_t, size_t> cpdim {p.size(), p.front().size()};

    CPointsType cp(cpdim.first * cpdim.second);
    for (size_t i = 0; i < cpdim.first; ++i) {
        for (size_t j = 0; j < cpdim.second; ++j) {
            cp[j + cpdim.second * i] = Super::CPoint(p[i][j].raw(), _w(i, j));
        }
    }

    c_ = Super(order, cpdim, k, std::move(cp));
}

void BSplineSurface::Impl::init_surface(
    std::pair<size_t, size_t> order,
    std::pair<std::vector<double>, std::vector<double>>&& k,
    const std::vector<std::vector<Point>>& p,
    const std::vector<std::vector<double>>& w)
{
    std::pair<size_t, size_t> cpdim {p.size(), p.front().size()};

    CPointsType cp(cpdim.first * cpdim.second);
    for (size_t i = 0; i < cpdim.first; ++i) {
        for (size_t j = 0; j < cpdim.second; ++j) {
            cp[j + cpdim.second * i] = Super::CPoint(p[i][j].raw(), _w(i, j));
        }
    }

    c_ = Super(order, cpdim, std::move(k), std::move(cp));
}

#undef _w

} // namespace gm