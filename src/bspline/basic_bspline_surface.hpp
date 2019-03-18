#ifndef GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_HPP_
#define GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_HPP_

#include <bspline/basic_bspline_curve.hpp>
#include <bspline/cox_de_boor.hpp>
#include <bspline/util.hpp>
#include <bspline/wpoint.hpp>
#include <gm/compare.hpp>
#include <gm/surf_point.hpp>

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
    enum class RefineDir { U_DIRECTION = 0, V_DIRECTION = 1 };

    BasicBSplineSurface()
        : order_(0, 0)
        , shape_(order_)
        , knots_()
        , cpoints_()
        , cdb_()
    {
    }

    BasicBSplineSurface(const BezierSurfacePatch<N>& patch)
        : order_(patch.order())
        , shape_(order_)
        , knots_(std::move(patch.knots()))
        , cpoints_(patch.cpoints())
        , cdb_()
    {
        init_cdb();
    }

    BasicBSplineSurface(BezierSurfacePatch<N>&& patch)
        : order_(patch.order())
        , shape_(order_)
        , knots_(std::move(patch.knots()))
        , cpoints_(std::move(patch.cpoints()))
        , cdb_()
    {
        init_cdb();
    }

    BasicBSplineSurface(std::pair<size_t, size_t> order,
                        std::pair<size_t, size_t> cpdim,
                        const KnotsType& knots, CPointsType&& cpoints)
        : order_(order)
        , shape_(cpdim)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_()
    {
        init_cdb();
    }

    BasicBSplineSurface(std::pair<size_t, size_t> order,
                        std::pair<size_t, size_t> cpdim, KnotsType&& knots,
                        CPointsType&& cpoints)
        : order_(order)
        , shape_(cpdim)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_()
    {
        init_cdb();
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

    const std::pair<size_t, size_t>& shape() const noexcept
    {
        return shape_;
    }

    const CPointsType& cpoints() const noexcept
    {
        return cpoints_;
    }

    const CPoint& operator[](const std::pair<size_t, size_t>& i) const noexcept
    {
        return cpoints_[ind(i.first, i.second)];
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

    std::vector<BezierSurfacePatch<N>> bezier_patches() const noexcept
    {
        throw std::runtime_error("not implemented");
        return {};
    }

    BasicBSplineSurface& refine_knots(const KnotsType& to_insert)
    {
        u_refine_knots(to_insert.first);
        v_refine_knots(to_insert.second);

        return *this;
    }

    BasicBSplineSurface& u_refine_knots(const std::vector<double>& to_insert)
    {
        static constexpr auto npos = size_t(-1);

        if (to_insert.empty()) {
            return *this;
        }

        auto& knots = knots_.first;
        auto& order = order_.first;
        auto& shape = shape_.first;
        auto& second = shape_.second;

        auto a = find_span(to_insert.front(), order, knots);
        auto b = find_span(to_insert.back(), order, knots);

        auto p = order - 1;
        auto n = shape - 1;
        auto m = knots.size() - 1;
        auto s = to_insert.size();
        auto i = b + p - 1;
        auto k = b + p + s - 1;
        size_t j, l, r;

        std::vector<CPoint> new_cpoints((shape + s) * second);
        std::vector<double> new_knots(knots.size() + s);
        for (j = 0; j <= a; ++j) {
            new_knots[j] = knots[j];
        }
        for (j = b + p; j <= m; ++j) {
            new_knots[j + s] = knots[j];
        }
        for (j = 0; j <= a - p; ++j) {
            for (r = 0; r < second; ++r) {
                new_cpoints[ind(j, r)] = cpoints_[ind(j, r)];
            }
        }
        for (j = b - 1; j <= n; ++j) {
            for (r = 0; r < second; ++r) {
                new_cpoints[ind(j + s, r)] = cpoints_[ind(j, r)];
            }
        }

        for (j = s - 1; j != npos; --j) {
            while (gm::cmp::le(to_insert[j], knots[i]) && i > a) {
                for (r = 0; r < second; ++r) {
                    new_cpoints[ind(k - p - 1, r)]
                        = cpoints_[ind(i - p - 1, r)];
                }
                new_knots[k] = knots[i];
                --k;
                --i;
            }
            for (r = 0; r < second; ++r) {
                new_cpoints[ind(k - p - 1, r)] = new_cpoints[ind(k - p, r)];
            }

            for (l = 1; l <= p; ++l) {
                auto pos = k - p + l;
                auto alpha = new_knots[k + l] - to_insert[j];

                if (gm::cmp::zero(alpha)) {
                    for (r = 0; r < second; ++r) {
                        auto& cur = new_cpoints[ind(pos, r)];
                        auto& prv = new_cpoints[ind(pos - 1, r)];

                        prv = cur;
                    }
                } else {
                    alpha /= new_knots[k + l] - knots[i - p + l];

                    for (r = 0; r < second; ++r) {
                        auto& cur = new_cpoints[ind(pos, r)];
                        auto& prv = new_cpoints[ind(pos - 1, r)];

                        prv = alpha * prv + (1. - alpha) * cur;
                    }
                }
            }
            new_knots[k] = to_insert[j];
            --k;
        }

        shape_ = std::make_pair(shape + s, second);
        knots_.first = std::move(new_knots);
        cpoints_ = std::move(new_cpoints);
        cdb_.clear();
        init_cdb();

        return *this;
    }

    BasicBSplineSurface& v_refine_knots(const std::vector<double>& to_insert)
    {
        static constexpr auto npos = size_t(-1);

        if (to_insert.empty()) {
            return *this;
        }

        auto& knots = knots_.second;
        auto& order = order_.second;
        auto& shape = shape_.second;
        auto& first = shape_.first;

        auto a = find_span(to_insert.front(), order, knots);
        auto b = find_span(to_insert.back(), order, knots);

        auto p = order - 1;
        auto n = shape - 1;
        auto m = knots.size() - 1;
        auto s = to_insert.size();
        auto i = b + p - 1;
        auto k = b + p + s - 1;
        size_t j, l, r;

        std::vector<CPoint> new_cpoints(first * (shape + s));
        std::vector<double> new_knots(knots.size() + s);
        for (j = 0; j <= a; ++j) {
            new_knots[j] = knots[j];
        }
        for (j = b + p; j <= m; ++j) {
            new_knots[j + s] = knots[j];
        }
        for (r = 0; r < first; ++r) {
            for (j = 0; j <= a - p; ++j) {
                new_cpoints[ind(r, j)] = cpoints_[ind(r, j)];
            }
        }
        for (r = 0; r < first; ++r) {
            for (j = b - 1; j <= n; ++j) {
                new_cpoints[ind(r, j + s)] = cpoints_[ind(r, j)];
            }
        }

        for (j = s - 1; j != npos; --j) {
            while (gm::cmp::le(to_insert[j], knots[i]) && i > a) {
                for (r = 0; r < first; ++r) {
                    new_cpoints[ind(r, k - p - 1)]
                        = cpoints_[ind(r, i - p - 1)];
                }
                new_knots[k] = knots[i];
                --k;
                --i;
            }
            for (r = 0; r < first; ++r) {
                new_cpoints[ind(r, k - p - 1)] = new_cpoints[ind(r, k - p)];
            }

            for (l = 1; l <= p; ++l) {
                auto pos = k - p + l;
                auto alpha = new_knots[k + l] - to_insert[j];

                if (gm::cmp::zero(alpha)) {
                    for (r = 0; r < first; ++r) {
                        auto& cur = new_cpoints[ind(r, pos)];
                        auto& prv = new_cpoints[ind(r, pos - 1)];

                        prv = cur;
                    }
                } else {
                    alpha /= new_knots[k + l] - knots[i - p + l];

                    for (r = 0; r < first; ++r) {
                        auto& cur = new_cpoints[ind(r, pos)];
                        auto& prv = new_cpoints[ind(r, pos - 1)];

                        prv = alpha * prv + (1. - alpha) * cur;
                    }
                }
            }
            new_knots[k] = to_insert[j];
            --k;
        }

        shape_ = std::make_pair(first, shape + s);
        knots_.second = std::move(new_knots);
        cpoints_ = std::move(new_cpoints);
        cdb_.clear();
        init_cdb();

        return *this;
    }

protected:
    size_t ind(size_t i, size_t j) const noexcept
    {
        return j + shape_.second * i;
    }

private:
    void init_cdb()
    {
        cdb_.reserve(shape_.first);
        VectorView<double> knots(knots_.second);

        auto ptr = cpoints_.data();
        auto end = &cpoints_.back();
        for (; ptr < end; ptr += shape_.second) {
            cdb_.emplace_back(order_.second, knots,
                              VectorView<CPoint>(ptr, shape_.second));
        }
        if (ptr == end) {
            end = nullptr;
        }
    }

    std::pair<size_t, size_t> order_;
    std::pair<size_t, size_t> shape_;
    KnotsType knots_;
    CPointsType cpoints_;
    CoxDeBoorType cdb_;
};

#endif // GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_HPP_