#include "cox_de_boor.h"

#include <util/debug.h>
#include <util/math.h>

using namespace std;

CoxDeBoor::CoxDeBoor()
    : order(0)
    , knots()
    , points()
    , weights()
{
}

CoxDeBoor::CoxDeBoor(size_t p_order, const vector<double>& p_knots,
                     const vector<Point>& p_points,
                     const vector<double>& p_weights)
    : order(p_order)
    , knots(p_knots)
    , points(p_points)
    , weights(p_weights)
{
}

CoxDeBoor::CoxDeBoorProxy CoxDeBoor::get_proxy(double u) const
{
    u = pad(u);
    return CoxDeBoorProxy(u, get_interval(u), *this);
}

double CoxDeBoor::pad(double u) const
{
    return ::pad(u, knots.front(), knots.back());
}

size_t CoxDeBoor::get_interval(double u) const
{
    size_t result = 0;
    size_t end = knots.size() - 1;
    for (size_t i = 0; i < end; ++i) {
        if (knots[i] <= u && u < knots[i + 1]) {
            result = i;
            break;
        }
    }
    return result;
}

CoxDeBoor::ResultType CoxDeBoor::eval(double u, size_t p, vector<double>& w,
                                      vector<Point>& wr) const
{
    auto n = wr.size();
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = n - 1; j >= i; --j) {
            size_t k = p - n + 1 + j;
            auto d1 = u - knots[k];
            auto d2 = knots[k + n - i] - u;
            wr[j] = (wr[j] * d1 + wr[j - 1] * d2) / (d1 + d2);
            w[j] = (w[j] * d1 + w[j - 1] * d2) / (d1 + d2);
        }
    }
    return {wr.back(), w.back()};
}

double CoxDeBoor::get_w(size_t i, size_t k) const
{
    double result;
    if (k == 0) {
        result = weights[i];
    } else {
        result = double(order - k) * (get_w(i, k - 1) - get_w(i - 1, k - 1))
            / (knots[i + order - 1] - knots[i]);
    }
    return result;
}

Point CoxDeBoor::get_wr(size_t i, size_t k) const
{
    Point result;
    if (k == 0) {
        result = points[i] * weights[i];
    } else {
        result = double(order - k) * (get_wr(i, k - 1) - get_wr(i - 1, k - 1))
            / (knots[i + order - 1] - knots[i]);
    }
    return result;
}

CoxDeBoor::CoxDeBoorProxy::CoxDeBoorProxy(double pu, size_t pp,
                                          const CoxDeBoor& parent)
    : u(pu)
    , p(pp)
    , parent_(parent)
{
}

CoxDeBoor::ResultType CoxDeBoor::CoxDeBoorProxy::get(size_t k) const
{
    ResultType result;

    if (k >= parent_.order) {
        result = {Point(), 0};
    } else {
        auto n = parent_.order - k;
        vector<double> w(n);
        vector<Point> wr(n);

        for (size_t i = 0; i < n; ++i) {
            w[i] = parent_.get_w(p - n + 1 + i, k);
            wr[i] = parent_.get_wr(p - n + 1 + i, k);
        }
        result = parent_.eval(u, p, w, wr);
    }

    return result;
}

vector<CoxDeBoor::ResultType> CoxDeBoor::CoxDeBoorProxy::range(size_t k) const
{
    vector<ResultType> result(k);
    for (size_t i = 0; i < k; ++i) {
        result[i] = get(i);
    }
    return result;
}

Point CoxDeBoor::CoxDeBoorProxy::f() const
{
    auto q = get(0);
    return q.r / q.w;
}

Vec CoxDeBoor::CoxDeBoorProxy::df() const
{
    auto q = range(2);
    return Vec(q[1].r / q[0].w - q[1].w * q[0].r / sqr(q[0].w));
}

Vec CoxDeBoor::CoxDeBoorProxy::df2() const
{
    auto q = range(3);
    return Vec(-2 * q[1].r * q[1].w / q[0].w
               + 2 * q[0].r * q[1].w / sqr(q[0].w) + q[2].r
               - q[2].w * q[0].r / q[0].w)
        / q[0].w;
}

pair<Point, double> cox_de_boor(double u, size_t p, size_t order,
                                const vector<double>& knots,
                                const vector<Point>& points,
                                const vector<double>& weights)
{
    vector<double> w(order);
    vector<Point> r(order);

    size_t i, j, k;
    double d1, d2;

    for (i = 0; i < order; ++i) {
        w[i] = weights[p - order + 1 + i];
        r[i] = points[p - order + 1 + i] * w[i];
    }
    for (i = 1; i < order; ++i) {
        for (j = order - 1; j >= i; --j) {
            k = p - order + 1 + j;
            d1 = u - knots[k];
            d2 = knots[k + order - i] - u;
            r[j] = (r[j] * d1 + r[j - 1] * d2) / (d1 + d2);
            w[j] = (w[j] * d1 + w[j - 1] * d2) / (d1 + d2);
        }
    }
    return {r.back() / w.back(), w.back()};
}

size_t get_interval(double u, const vector<double>& knots)
{
    size_t end = knots.size() - 1;
    for (size_t i = 0; i != end; ++i)
        if (knots[i] <= u && u < knots[i + 1])
            return i;
    return 0;
}