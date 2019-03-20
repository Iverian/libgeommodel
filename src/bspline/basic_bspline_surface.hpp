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
#include <utility>
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

    BezierSurfacePatch(std::pair<size_t, size_t> order,
                       std::pair<gm::SurfPoint, gm::SurfPoint> limits,
                       CPointsType&& cpoints)
        : order_(std::move(order))
        , limits_(std::move(limits))
        , cpoints_(cpoints)
    {
    }

    const std::pair<size_t, size_t>& order() const noexcept
    {
        return order_;
    }
    gm::SurfPoint& pfront() noexcept
    {
        return limits_.first;
    }
    gm::SurfPoint& pback() noexcept
    {
        return limits_.second;
    }
    const gm::SurfPoint& pfront() const noexcept
    {
        return limits_.first;
    }
    const gm::SurfPoint& pback() const noexcept
    {
        return limits_.second;
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
    std::pair<gm::SurfPoint, gm::SurfPoint> limits_;
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

    struct PartialBezierPatches {
        std::vector<CPointsType> cpoints;
        std::vector<double> unique_knots;
        std::pair<size_t, size_t> shape;

        auto size() const noexcept
        {
            return cpoints.size();
        }
    };

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

    auto f(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).f().p();
    }

    auto dfu(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).df().p();
    }

    auto dfv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df(); });

        return ucdb.proxy(p.u).f().p();
    }

    auto dfuu(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).f(); });

        return ucdb.proxy(p.u).df2().p();
    }

    auto dfuv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df(); });

        return ucdb.proxy(p.u).df().p();
    }

    auto dfvv(const gm::SurfPoint& p) const noexcept
    {
        CPointsType ubuf(cdb_.size());
        CoxDeBoorTypeU ucdb(order_.first, knots_.first, ubuf);
        std::transform(std::begin(cdb_), std::end(cdb_), std::begin(ubuf),
                       [&p](const auto& x) { return x.proxy(p.v).df2(); });

        return ucdb.proxy(p.u).f().p();
    }

    auto bezier_patches() const
    {
        std::vector<BezierSurfacePatch<N>> result;
        auto u_patches = u_bezier_patches(cpoints_, shape_);
        auto n = u_patches.size();

        for (size_t i = 0; i < n; ++i) {
            auto v_patches
                = v_bezier_patches(u_patches.cpoints[i], u_patches.shape);
            auto m = v_patches.size();
            auto& ua = u_patches.unique_knots[i];
            auto& ub = u_patches.unique_knots[i + 1];

            for (size_t j = 0; j < m; ++j) {
                auto& cp = v_patches.cpoints[j];
                auto& va = v_patches.unique_knots[j];
                auto& vb = v_patches.unique_knots[j + 1];

                result.emplace_back(v_patches.shape,
                                    std::make_pair(gm::SurfPoint(ua, va),
                                                   gm::SurfPoint(ub, vb)),
                                    std::move(cp));
            }
        }

        return result;
    }

    auto refine_knots(const KnotsType& to_insert)
    {
        u_refine_knots(to_insert.first);
        v_refine_knots(to_insert.second);

        return *this;
    }

    auto u_refine_knots(const std::vector<double>& to_insert)
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
        auto new_shape = shape + s;
        size_t j, l, r;

        std::vector<CPoint> new_cpoints(new_shape * second);
        std::vector<double> new_knots(knots.size() + s);
        for (j = 0; j <= a; ++j) {
            new_knots[j] = knots[j];
        }
        for (j = b + p; j <= m; ++j) {
            new_knots[j + s] = knots[j];
        }
        for (j = 0; j <= a - p; ++j) {
            for (r = 0; r < second; ++r) {
                new_cpoints[second * j + r] = cpoints_[ind(j, r)];
            }
        }
        for (j = b - 1; j <= n; ++j) {
            for (r = 0; r < second; ++r) {
                new_cpoints[second * (j + s) + r] = cpoints_[ind(j, r)];
            }
        }

        for (j = s - 1; j != npos; --j) {
            while (gm::cmp::le(to_insert[j], knots[i]) && i > a) {
                for (r = 0; r < second; ++r) {
                    new_cpoints[second * (k - p - 1) + r]
                        = cpoints_[ind(i - p - 1, r)];
                }
                new_knots[k] = knots[i];
                --k;
                --i;
            }
            for (r = 0; r < second; ++r) {
                new_cpoints[second * (k - p - 1) + r]
                    = new_cpoints[second * (k - p) + r];
            }

            for (l = 1; l <= p; ++l) {
                auto pos = k - p + l;
                auto alpha = new_knots[k + l] - to_insert[j];

                if (gm::cmp::zero(alpha)) {
                    for (r = 0; r < second; ++r) {
                        auto& cur = new_cpoints[second * pos + r];
                        auto& prv = new_cpoints[second * (pos - 1) + r];

                        prv = cur;
                    }
                } else {
                    alpha /= new_knots[k + l] - knots[i - p + l];

                    for (r = 0; r < second; ++r) {
                        auto& cur = new_cpoints[second * pos + r];
                        auto& prv = new_cpoints[second * (pos - 1) + r];

                        prv = alpha * prv + (1. - alpha) * cur;
                    }
                }
            }
            new_knots[k] = to_insert[j];
            --k;
        }

        shape_ = std::make_pair(new_shape, second);
        knots_.first = std::move(new_knots);
        cpoints_ = std::move(new_cpoints);
        cdb_.clear();
        init_cdb();

        return *this;
    }

    auto v_refine_knots(const std::vector<double>& to_insert)
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
        auto new_shape = shape + s;
        size_t j, l, r;

        std::vector<CPoint> new_cpoints(first * new_shape);
        std::vector<double> new_knots(knots.size() + s);
        for (j = 0; j <= a; ++j) {
            new_knots[j] = knots[j];
        }
        for (j = b + p; j <= m; ++j) {
            new_knots[j + s] = knots[j];
        }
        for (r = 0; r < first; ++r) {
            for (j = 0; j <= a - p; ++j) {
                new_cpoints[new_shape * r + j] = cpoints_[ind(r, j)];
            }
        }
        for (r = 0; r < first; ++r) {
            for (j = b - 1; j <= n; ++j) {
                new_cpoints[new_shape * r + j + s] = cpoints_[ind(r, j)];
            }
        }

        for (j = s - 1; j != npos; --j) {
            while (gm::cmp::le(to_insert[j], knots[i]) && i > a) {
                for (r = 0; r < first; ++r) {
                    new_cpoints[new_shape * r + k - p - 1]
                        = cpoints_[ind(r, i - p - 1)];
                }
                new_knots[k] = knots[i];
                --k;
                --i;
            }
            for (r = 0; r < first; ++r) {
                new_cpoints[new_shape * r + k - p - 1]
                    = new_cpoints[new_shape * r + k - p];
            }

            for (l = 1; l <= p; ++l) {
                auto pos = k - p + l;
                auto alpha = new_knots[k + l] - to_insert[j];

                if (gm::cmp::zero(alpha)) {
                    for (r = 0; r < first; ++r) {
                        auto& cur = new_cpoints[new_shape * r + pos];
                        auto& prv = new_cpoints[new_shape * r + pos - 1];

                        prv = cur;
                    }
                } else {
                    alpha /= new_knots[k + l] - knots[i - p + l];

                    for (r = 0; r < first; ++r) {
                        auto& cur = new_cpoints[new_shape * r + pos];
                        auto& prv = new_cpoints[new_shape * r + pos - 1];

                        prv = alpha * prv + (1. - alpha) * cur;
                    }
                }
            }
            new_knots[k] = to_insert[j];
            --k;
        }

        shape_ = std::make_pair(first, new_shape);
        knots_.second = std::move(new_knots);
        cpoints_ = std::move(new_cpoints);
        cdb_.clear();
        init_cdb();

        return *this;
    }

protected:
    auto ind(size_t i, size_t j) const noexcept
    {
        return j + shape_.second * i;
    }

private:
    auto u_bezier_patches(const CPointsType& cp,
                          std::pair<size_t, size_t> shape) const
    {
        static constexpr auto npos = size_t(-1);

        PartialBezierPatches result;
        std::vector<double> alphas(order_.first);

        auto& knots = knots_.first;
        auto p = order_.first - 1;
        auto m = knots.size() - 1;
        auto size = order_.first * shape.second;
        auto a = p;
        auto b = p + 1;
        size_t i, j, k, l, nb = 0;

        result.shape = std::make_pair(order_.first, shape.second);
        result.cpoints.emplace_back(size);
        for (i = 0; i < order_.first; ++i) {
            for (j = 0; j < shape.second; ++j) {
                result.cpoints[nb][shape.second * i + j]
                    = cp[shape.second * i + j];
            }
        }

        result.unique_knots.emplace_back(knots[a]);
        while (b < m) {
            i = b;
            for (; b < m && gm::cmp::near(knots[b], knots[b + 1]); ++b)
                ;
            result.unique_knots.emplace_back(knots[b]);

            if (b < m) {
                result.cpoints.emplace_back(size);
            }

            auto mult = b - i + 1;
            if (mult < p) {
                auto numer = knots[b] - knots[a];
                auto r = p - mult;
                for (j = p; j > mult; --j) {
                    alphas[j - mult - 1] = numer / (knots[a + j] - knots[a]);
                }
                for (j = 1; j < r + 1; ++j) {
                    auto save = r - j;
                    auto s = mult + j;
                    for (k = p; k >= s && k != npos; --k) {
                        auto& alpha = alphas[k - s];

                        for (l = 0; l < shape.second; ++l) {
                            auto& cur
                                = result.cpoints[nb][shape.second * k + l];
                            auto& prv
                                = result
                                      .cpoints[nb][shape.second * (k - 1) + l];

                            cur = alpha * cur + (1. - alpha) * prv;
                        }
                    }

                    if (b < m) {
                        for (l = 0; l < shape.second; ++l) {
                            result.cpoints[nb + 1][shape.second * save + l]
                                = result.cpoints[nb][shape.second * p + l];
                        }
                    }
                }
                std::fill(std::begin(alphas), std::end(alphas), 0.);
            }

            ++nb;
            if (b < m) {
                for (i = p - mult; i < order_.first; ++i) {
                    for (j = 0; j < shape.second; ++j) {
                        result.cpoints[nb][shape.second * i + j]
                            = cp[shape.second * (b - p + i) + j];
                    }
                }
                a = b;
                ++b;
            }
        }

        return result;
    }

    auto v_bezier_patches(const CPointsType& cp,
                          std::pair<size_t, size_t> shape) const
    {
        static constexpr auto npos = size_t(-1);

        PartialBezierPatches result;
        std::vector<double> alphas(order_.second);

        auto& knots = knots_.second;
        auto p = order_.second - 1;
        auto m = knots.size() - 1;
        auto size = shape.first * order_.second;
        auto a = p;
        auto b = p + 1;
        size_t i, j, k, l, nb = 0;

        result.shape = std::make_pair(shape.first, order_.second);
        result.cpoints.emplace_back(size);
        for (i = 0; i < shape.first; ++i) {
            for (j = 0; j < order_.second; ++j) {
                result.cpoints[nb][order_.second * i + j]
                    = cp[shape.second * i + j];
            }
        }

        result.unique_knots.emplace_back(knots[a]);
        while (b < m) {
            i = b;
            for (; b < m && gm::cmp::near(knots[b], knots[b + 1]); ++b)
                ;
            result.unique_knots.emplace_back(knots[b]);

            if (b < m) {
                result.cpoints.emplace_back(size);
            }

            auto mult = b - i + 1;
            if (mult < p) {
                auto numer = knots[b] - knots[a];
                auto r = p - mult;
                for (j = p; j > mult; --j) {
                    alphas[j - mult - 1] = numer / (knots[a + j] - knots[a]);
                }
                for (j = 1; j < r + 1; ++j) {
                    auto save = r - j;
                    auto s = mult + j;
                    for (k = p; k >= s && k != npos; --k) {
                        auto& alpha = alphas[k - s];

                        for (l = 0; l < shape.first; ++l) {
                            auto& cur
                                = result.cpoints[nb][order_.second * l + k];
                            auto& prv = result.cpoints[nb][order_.second * l
                                                           + (k - 1)];

                            cur = alpha * cur + (1. - alpha) * prv;
                        }
                    }

                    if (b < m) {
                        for (l = 0; l < shape.first; ++l) {
                            result.cpoints[nb + 1][order_.second * l + save]
                                = result.cpoints[nb][order_.second * l + p];
                        }
                    }
                }
                std::fill(std::begin(alphas), std::end(alphas), 0.);
            }

            ++nb;
            if (b < m) {
                for (i = 0; i < shape.first; ++i) {
                    for (j = p - mult; j < order_.second; ++j) {
                        result.cpoints[nb][order_.second * i + j]
                            = cp[shape.second * i + (b - p + j)];
                    }
                }
                a = b;
                ++b;
            }
        }

        return result;
    }

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
    }

    std::pair<size_t, size_t> order_;
    std::pair<size_t, size_t> shape_;
    KnotsType knots_;
    CPointsType cpoints_;
    CoxDeBoorType cdb_;
};

#endif // GEOMMODEL_SRC_BASIC_BSPLINE_SURFACE_HPP_