#ifndef GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_H_
#define GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_H_

#include <bspline/basic_bspline_curve.h>
#include <bspline/cox_de_boor.h>
#include <bspline/wpoint.h>
#include <gm/compare.h>
#include <gm/surf_point.h>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <vector>

template <size_t N>
struct BSplineSurfaceTraits {
    using CPoint = gm::WPoint<double, N>;
    using CoxDeBoorTypeU = gm::CoxDeBoor<double, N>;
    using CoxDeBoorType = std::vector<CoxDeBoorTypeU>;
    using KnotsType = std::pair<std::vector<double>, std::vector<double>>;
    using CPointsType = std::vector<CPoint>;
};

template <size_t N>
class BezierSurfacePatch {
public:
    using CPoint = typename BSplineSurfaceTraits<N>::CPoint;
    using KnotsType = typename BSplineSurfaceTraits<N>::KnotsType;
    using CPointsType = typename BSplineSurfaceTraits<N>::CPointsType;

    explicit BezierSurfacePatch(std::pair<size_t, size_t> order = {0, 0})
        : order_(order)
        , limits_({{0, 0}, {1, 1}})
        , cpoints_(order.first * order.second)
    {
    }

    const std::pair<size_t, size_t>& order() const noexcept
    {
        return order_;
    }
    gm::SurfPoint& pfront() noexcept
    {
        return limits_[0];
    }
    gm::SurfPoint& pback() noexcept
    {
        return limits_[1];
    }
    const gm::SurfPoint& pfront() const noexcept
    {
        return limits_[0];
    }
    const gm::SurfPoint& pback() const noexcept
    {
        return limits_[1];
    }
    KnotsType knots() const noexcept
    {
        KnotsType result;
        auto& [su, sv] = order_;

        result.first.resize(2 * su);
        for (size_t i = 0; i < su; ++i) {
            result.first[i] = pfront().u;
            result.first[su + i] = pback().u;
        }

        result.second.resize(2 * sv);
        for (size_t i = 0; i < sv; ++i) {
            result.second[i] = pfront().v;
            result.second[sv + i] = pback().v;
        }

        return result;
    }
    CPointsType& cpoints() noexcept
    {
        return cpoints_;
    }
    const CPointsType& cpoints() const noexcept
    {
        return cpoints_;
    }

    CPoint& operator[](const std::pair<size_t, size_t>& i) noexcept
    {
        return cpoints_[i.second + order_.second * i.first];
    }

    const CPoint& operator[](const std::pair<size_t, size_t>& i) const noexcept
    {
        return cpoints_[i.second + order_.second * i.first];
    }

private:
    std::pair<size_t, size_t> order_;
    std::array<gm::SurfPoint, 2> limits_;
    CPointsType cpoints_;
};

template <size_t N>
class BasicBSplineSurface {
    using Traits = BSplineSurfaceTraits<N>;

public:
    using CPoint = typename Traits::CPoint;
    using KnotsType = typename Traits::KnotsType;
    using CPointsType = typename Traits::CPointsType;
    using CoxDeBoorType = typename Traits::CoxDeBoorType;
    using CoxDeBoorTypeU = typename Traits::CoxDeBoorTypeU;

    BasicBSplineSurface()
        : order_(0, 0)
        , knots_()
        , cpoints_()
        , cdb_()
    {
    }

    BasicBSplineSurface(const BezierSurfacePatch<N>& patch)
        : order_(patch.order())
        , knots_(std::move(patch.knots()))
        , cpoints_(patch.cpoints())
        , cdb_()
    {
        init_cdb(order_);
    }

    BasicBSplineSurface(BezierSurfacePatch<N>&& patch)
        : order_(patch.order())
        , knots_(std::move(patch.knots()))
        , cpoints_(std::move(patch.cpoints()))
        , cdb_()
    {
        init_cdb(order_);
    }

    BasicBSplineSurface(std::pair<size_t, size_t> order,
                        std::pair<size_t, size_t> cpdim,
                        const KnotsType& knots, CPointsType&& cpoints)
        : order_(order)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_()
    {
        init_cdb(cpdim);
    }

    BasicBSplineSurface(std::pair<size_t, size_t> order,
                        std::pair<size_t, size_t> cpdim, KnotsType&& knots,
                        CPointsType&& cpoints)
        : order_(order)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_()
    {
        init_cdb(cpdim);
    }

    const std::pair<size_t, size_t>& order() const noexcept
    {
        return order_;
    }

    gm::SurfPoint pfront() const noexcept
    {
        return {knots_.first.front(), knots_.second.front()};
    }

    gm::SurfPoint pback() const noexcept
    {
        return {knots_.first.back(), knots_.second.back()};
    }

    const KnotsType& knots() const noexcept
    {
        return knots_;
    }

    std::pair<size_t, size_t> shape() const noexcept
    {
        auto n = cdb_.size();
        auto m = cpoints_.size() / n;

        return {n, m};
    }

    const CPointsType& cpoints() const noexcept
    {
        return cpoints_;
    }

    const CPoint& operator[](const std::pair<size_t, size_t>& i) const noexcept
    {
        auto m = cpoints_.size() / cdb_.size();
        return cpoints_[i.second + m * i.first];
    }

    typename CPoint::Proj f(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).f().p();
    }

    typename CPoint::Proj dfu(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).df().p();
    }

    typename CPoint::Proj dfv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df(); });

        return ucdb.proxy(p.u).f().p();
    }

    typename CPoint::Proj dfuu(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).df2().p();
    }

    typename CPoint::Proj dfuv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df(); });

        return ucdb.proxy(p.u).df().p();
    }

    typename CPoint::Proj dfvv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df2(); });

        return ucdb.proxy(p.u).f().p();
    }

private:
    void init_cdb(const std::pair<size_t, size_t>& cpdim) noexcept
    {
        cdb_.reserve(cpdim.first);
        VectorView<double> knots(knots_.second);

        auto ptr = cpoints_.data();
        auto end = &cpoints_.back();
        for (; ptr < end; ptr += cpdim.second) {
            cdb_.emplace_back(order_.second, knots,
                              VectorView<CPoint>(ptr, cpdim.second));
        }
        if (ptr == end) {
            end = nullptr;
        }
    }

    std::pair<size_t, size_t> order_;
    KnotsType knots_;
    CPointsType cpoints_;
    CoxDeBoorType cdb_;
};

#endif // GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_H_