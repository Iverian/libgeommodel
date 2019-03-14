#ifndef GEOMMODEL_SRC_BASIC_BSPLINE_CURVE_H_
#define GEOMMODEL_SRC_BASIC_BSPLINE_CURVE_H_

#include <bspline/cox_de_boor.h>
#include <bspline/wpoint.h>
#include <gm/compare.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

template <size_t N>
struct BSplineCurveTraits {
    using CPoint = gm::WPoint<double, N>;
    using CoxDeBoorType = gm::CoxDeBoor<double, N>;
    using KnotsType = std::vector<double>;
    using CPointsType = std::vector<CPoint>;
};

template <size_t N>
class BezierPatch {
public:
    using CPoint = typename BSplineCurveTraits<N>::CPoint;
    using KnotsType = typename BSplineCurveTraits<N>::KnotsType;
    using CPointsType = typename BSplineCurveTraits<N>::CPointsType;

    explicit BezierPatch(size_t order = 0)
        : pfront_(0)
        , pback_(1)
        , cpoints_(order)
    {
    }

    size_t order() const noexcept
    {
        return cpoints_.size();
    }
    double& pfront() noexcept
    {
        return pfront_;
    }
    double& pback() noexcept
    {
        return pback_;
    }
    CPointsType& cpoints() noexcept
    {
        return cpoints_;
    }
    const double& pfront() const noexcept
    {
        return pfront_;
    }
    const double& pback() const noexcept
    {
        return pback_;
    }
    const CPointsType& cpoints() const noexcept
    {
        return cpoints_;
    }

    KnotsType knots() const noexcept
    {
        auto s = order();
        KnotsType result(2 * s);
        for (decltype(s) i = 0; i < s; ++i) {
            result[i] = pfront_;
            result[s + i] = pback_;
        }
        return result;
    }

private:
    double pfront_;
    double pback_;
    CPointsType cpoints_;
};

template <size_t N>
class BasicBsplineCurve {
public:
    using CPoint = typename BSplineCurveTraits<N>::CPoint;
    using KnotsType = typename BSplineCurveTraits<N>::KnotsType;
    using CPointsType = typename BSplineCurveTraits<N>::CPointsType;
    using CoxDeBoorType = typename BSplineCurveTraits<N>::CoxDeBoorType;

    BasicBsplineCurve()
        : order_(0)
        , knots_()
        , cpoints_()
        , cdb_()
    {
    }

    BasicBsplineCurve(const BezierPatch<N>& patch)
        : order_(patch.order())
        , knots_(std::move(patch.knots()))
        , cpoints_(patch.cpoints())
        , cdb_(order_, knots_, cpoints_)
    {
    }

    BasicBsplineCurve(BezierPatch<N>&& patch)
        : order_(patch.order())
        , knots_(std::move(patch.knots()))
        , cpoints_(std::move(patch.cpoints()))
        , cdb_(order_, knots_, cpoints_)
    {
    }

    BasicBsplineCurve(size_t order, KnotsType&& knots, CPointsType&& cpoints)
        : order_(order)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_(order_, knots_, cpoints_)
    {
    }

    BasicBsplineCurve(size_t order, const KnotsType& knots,
                      const CPointsType& cpoints)
        : order_(order)
        , knots_(knots)
        , cpoints_(cpoints)
        , cdb_(order_, knots_, cpoints_)
    {
    }

    typename CPoint::Proj f(double u) const noexcept
    {
        return cdb_.proxy(u).f().p();
    }

    typename CPoint::Proj df(double u) const noexcept
    {
        return cdb_.proxy(u).df().p();
    }

    typename CPoint::Proj df2(double u) const noexcept
    {
        return cdb_.proxy(u).df2().p();
    }

    double pfront() const noexcept
    {
        return knots_.front();
    }
    double pback() const noexcept
    {
        return knots_.back();
    }

    size_t order() const noexcept
    {
        return order_;
    }
    const KnotsType& knots() const noexcept
    {
        return knots_;
    }
    const CPointsType& cpoints() const noexcept
    {
        return cpoints_;
    }

    BasicBsplineCurve refine_knots(const KnotsType& to_insert)
    {
        static constexpr auto npos = size_t(-1);

        auto p = order_ - 1;
        auto n = cpoints_.size() - 1;
        auto m = knots_.size() - 1;
        auto s = to_insert.size();
        auto r = s - 1;

        auto a = cdb_.interval(to_insert.front());
        auto b = cdb_.interval(to_insert.back()) + 1;
        size_t j, l;

        std::vector<CPoint> new_cpoints(cpoints_.size() + s);
        std::vector<double> new_knots(knots_.size() + s);

        for (j = 0; j <= a - p; ++j) {
            new_cpoints[j] = cpoints_[j];
        }
        for (j = b - 1; j <= n; ++j) {
            new_cpoints[j + r + 1] = cpoints_[j];
        }
        for (j = 0; j <= a; ++j) {
            new_knots[j] = knots_[j];
        }
        for (j = b + p; j <= m; ++j) {
            new_knots[j + r + 1] = knots_[j];
        }

        auto i = b + p - 1;
        auto k = b + p + r;
        for (j = r; j != npos; --j) {
            while (gm::cmp::le(to_insert[j], knots_[i]) && i > a) {
                new_cpoints[k - p - 1] = cpoints_[i - p - 1];
                new_knots[k] = knots_[i];
                --k;
                --i;
            }
            new_cpoints[k - p - 1] = new_cpoints[k - p];

            for (l = 1; l <= p; ++l) {
                auto ind = k - p + l;
                auto alpha = new_knots[k + l] - to_insert[j];
                auto& cur = new_cpoints[ind];
                auto& prev = new_cpoints[ind - 1];

                if (gm::cmp::zero(alpha)) {
                    prev = cur;
                } else {
                    alpha /= new_knots[k + l] - knots_[i - p + l];
                    prev = alpha * prev + (1. - alpha) * cur;
                }
            }
            new_knots[k] = to_insert[j];
            --k;
        }

        return BasicBsplineCurve<N>(order_, std::move(new_knots),
                                    std::move(new_cpoints));
    }
    std::vector<BezierPatch<N>> bezier_patches() const
    {
        static constexpr auto npos = size_t(-1);

        std::vector<BezierPatch<N>> result;
        std::vector<double> alphas(order_);

        auto p = order_ - 1;
        auto m = knots_.size() - 1;
        auto a = p;
        auto b = p + 1;
        size_t i, j, k, nb = 0;

        result.emplace_back(order_);
        for (i = 0; i < order_; ++i) {
            result[0].cpoints()[i] = cpoints_[i];
        }

        while (b < m) {
            i = b;
            for (; b < m && gm::cmp::near(knots_[b], knots_[b + 1]); ++b)
                ;
            result[nb].pfront() = knots_[a];
            result[nb].pback() = knots_[b];

            if (b < m) {
                result.emplace_back(order_);
            }

            auto mult = b - i + 1;
            if (mult < p) {
                auto numer = knots_[b] - knots_[a];
                auto r = p - mult;
                for (j = p; j > mult; --j) {
                    alphas[j - mult - 1] = numer / (knots_[a + j] - knots_[a]);
                }
                for (j = 1; j < r + 1; ++j) {
                    auto save = r - j;
                    auto s = mult + j;
                    for (k = p; k >= s && k != npos; --k) {
                        auto& alpha = alphas[k - s];
                        auto& cur = result[nb].cpoints()[k];
                        auto& prev = result[nb].cpoints()[k - 1];

                        cur = alpha * cur + (1 - alpha) * prev;
                    }

                    if (b < m) {
                        result[nb + 1].cpoints()[save]
                            = result[nb].cpoints()[p];
                    }
                }
                std::fill(std::begin(alphas), std::end(alphas), 0.);
            }

            ++nb;
            if (b < m) {
                for (i = p - mult; i < order_; ++i) {
                    result[nb].cpoints()[i] = cpoints_[b - p + i];
                }
                a = b;
                ++b;
            }
        }

        return result;
    }

private:
    size_t order_;
    KnotsType knots_;
    CPointsType cpoints_;
    CoxDeBoorType cdb_;
};

#endif // GEOMMODEL_SRC_BASIC_BSPLINE_CURVE_H_