#ifndef GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_HPP_
#define GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_HPP_

#include <bspline/wpoint.hpp>
#include <gm/compare.hpp>
#include <util/debug.hpp>
#include <util/vector_view.hpp>

#include <array>
#include <type_traits>
#include <vector>

namespace gm {

template <class T, size_t N>
class CoxDeBoor {
public:
    using scalar_type = std::decay_t<T>;
    using value_type = WPoint<scalar_type, N>;
    using reference = value_type&;
    using const_reference = const value_type&;

    friend class Proxy;

    class Proxy {
        friend class CoxDeBoor;

    public:
        Proxy(const CoxDeBoor& parent, scalar_type t)
            : parent_(&parent)
            , t_(t)
            , p_(parent_->interval(t_))
        {
        }

        value_type get(size_t k) const
        {
            value_type result;
            if (k >= parent_->order_) {
                result = value_type();
            } else {
                auto n = parent_->order_ - k;
                std::vector<value_type> cp(n);

                for (size_t i = 0; i < n; ++i) {
                    cp[i] = parent_->pget(p_ - n + 1 + i, k);
                }
                result = parent_->eval(t_, p_, cp);
            }

            return result;
        }

        std::vector<value_type> range(size_t k) const
        {
            std::vector<value_type> result(k);
            for (size_t i = 0; i < k; ++i) {
                result[i] = get(i);
            }
            return result;
        }

        value_type f() const
        {
            return get(0);
        }

        value_type df() const
        {
            auto r = range(2);
            auto p = r[1].wdiv(r[0].w()) - r[0].wmul(r[1].w() / sqr(r[0].w()));
            // auto p = r[1].wp() / r[0].w()
            //     - (r[1].w() / sqr(r[0].w())) * r[0].wp();

            return value_type(p.wp(), 1);
        }

        value_type df2() const
        {
            auto r = range(3);
            auto sw = sqr(r[0].w());
            auto p = r[2].wdiv(r[0].w()) - r[1].wmul(2 * r[1].w() / sw)
                + r[0].wmul((2 * r[1].w() / r[0].w() - r[2].w()) / sw);
            // auto p = (r[1].wp() * (-2 * r[1].w() / r[0].w())
            //           + r[0].wp() * (2 * r[1].w() / sqr(r[0].w())) +
            //           r[2].wp()
            //           - r[0].wp() * (r[2].w() / r[0].w()))
            //     / r[0].w();

            return value_type(p.wp(), 1);
        }

    private:
        const CoxDeBoor* parent_;
        scalar_type t_;
        size_t p_;
    };

    CoxDeBoor()
        : order_(0)
        , knots_()
        , cpoints_()
    {
    }

    CoxDeBoor(size_t order, const std::vector<scalar_type>& knots,
              const std::vector<value_type>& cpoints)
        : order_(order)
        , knots_(knots)
        , cpoints_(cpoints)
    {
    }

    CoxDeBoor(size_t order, VectorView<scalar_type> knots,
              VectorView<value_type> cpoints)
        : order_(std::move(order))
        , knots_(std::move(knots))
        , cpoints_(std::move(cpoints))
    {
    }

    Proxy proxy(scalar_type t) const
    {
        check_if(cmp::ge(t, knots_.front()) && gm::cmp::le(t, knots_.back()),
                 "Argument t = {} out of range [{}; {}]", t, knots_.front(),
                 knots_.back());

        return Proxy(*this, t);
    }

    size_t interval(scalar_type t) const
    {
        auto p = order_ - 1;
        auto n = knots_.size() - order_ - 1;

        if (gm::cmp::near(t, knots_[n + 1])) {
            return n;
        }

        auto low = p;
        auto high = n + 1;
        auto mid = (low + high) / 2;
        while (t < knots_[mid] || cmp::ge(t, knots_[mid + 1])) {
            if (t < knots_[mid]) {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low + high) / 2;
        }

        return mid;
    }

    value_type pget(size_t i, size_t k) const
    {
        value_type result;
        if (k == 0) {
            result = cpoints_.at(i);
        } else {
            result = scalar_type(order_ - k)
                * (pget(i, k - 1) - pget(i - 1, k - 1))
                / (knots_.at(i + order_ - 1) - knots_.at(i));
        }
        return result;
    }

    value_type eval(scalar_type t, size_t p, std::vector<value_type>& cp) const
    {
        auto n = cp.size();
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = n - 1; j >= i; --j) {
                auto k = p - n + 1 + j;
                auto a = t - knots_[k];
                auto b = knots_[k + n - i] - t;
                cp[j] = (cp[j] * a + cp[j - 1] * b) / (a + b);
            }
        }
        return cp.back();
    }

private:
    size_t order_;
    VectorView<scalar_type> knots_;
    VectorView<value_type> cpoints_;
};

} // namespace gm

#endif // GEOM_MODEL_SRC_GEOM_COX_DE_BOOR_HPP_