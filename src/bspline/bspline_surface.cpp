#include <bspline/bspline_surface_impl.hpp>
#include <gm/bspline_surface.hpp>

using namespace std;

namespace gm {

BSplineSurface::BSplineSurface()
    : pimpl_(make_unique<Impl>())
{
}

BSplineSurface::BSplineSurface(size_t du, size_t dv, const vector<double>& ku,
                               const vector<double>& kv,
                               const vector<vector<Point>>& p,
                               const vector<vector<double>>& w)
    : pimpl_(make_unique<Impl>(du, dv, ku, kv, p, w))
{
}

BSplineSurface::BSplineSurface(size_t du, size_t dv,
                               const vector<size_t>& ku_mult,
                               const vector<double>& ku_vals,
                               const vector<size_t>& kv_mult,
                               const vector<double>& kv_vals,
                               const vector<vector<Point>>& p,
                               const vector<vector<double>>& w)
    : pimpl_(
          make_unique<Impl>(du, dv, ku_mult, ku_vals, kv_mult, kv_vals, p, w))
{
}

Point BSplineSurface::f(const SurfPoint& p) const noexcept
{
    return pimpl_->f(p);
}

Vec BSplineSurface::dfu(const SurfPoint& p) const noexcept
{
    return pimpl_->dfu(p);
}

Vec BSplineSurface::dfv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfv(p);
}

Vec BSplineSurface::dfuu(const SurfPoint& p) const noexcept
{
    return pimpl_->dfuu(p);
}

Vec BSplineSurface::dfuv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfuv(p);
}

Vec BSplineSurface::dfvv(const SurfPoint& p) const noexcept
{
    return pimpl_->dfvv(p);
}

ostream& BSplineSurface::print(ostream& os) const
{
    return pimpl_->print(os);
}

SurfPoint BSplineSurface::project(const Point& p) const
{
    return pimpl_->project(p);
}

} // namespace gm

// struct BSplineSurface::Impl {
//     Impl(size_t du, size_t dv, const vector<double>& ku,
//          const vector<double>& kv, const vector<vector<Point>>& p,
//          const vector<vector<double>>& w = vector<vector<double>>());

//     Impl(size_t du, size_t dv, const vector<size_t>& ku_mult,
//          const vector<double>& ku_vals, const vector<size_t>& kv_mult,
//          const vector<double>& kv_vals, const vector<vector<Point>>& p,
//          const vector<vector<double>>& w = vector<vector<double>>());

//     Point f(const SurfPoint& p) const;
//     Vec dfu(const SurfPoint& p) const;
//     Vec dfv(const SurfPoint& p) const;
//     Vec dfuu(const SurfPoint& p) const;
//     Vec dfvv(const SurfPoint& p) const;
//     Vec dfuv(const SurfPoint& p) const;

//     ostream& print(ostream& os) const;

//     SurfPoint project(const Point& p, const BSplineSurface& c) const;

//     // pair<size_t, bool> is_init(const Point& p,
//     //                            const array<SurfPoint, 4>& control) const;

//     // optional<SurfPoint> project_u(const Point& p, double v, double
//     ubegin,
//     //                               double uend, double tol1, double tol2,
//     //                               size_t max_iter, size_t dcoeff,
//     //                               const BSplineSurface& c) const;
//     // optional<SurfPoint> project_v(const Point& p, double u, double
//     vbegin,
//     //                               double vend, double tol1, double tol2,
//     //                               size_t max_iter, size_t dcoeff,
//     //                               const BSplineSurface& c) const;

//     // SurfPoint project_s(const Point& p, const BSplineSurface& c) const;
//     // double ustep(const SurfPoint& p, const BSplineSurface& c) const;
//     // double vstep(const SurfPoint& p, const BSplineSurface& c) const;
//     // double init_between(const Point& p, const SurfPoint& a,
//     //                     const SurfPoint& b) const;

//     // optional<SurfPoint> newton_iter(const Point& p, const SurfPoint&
//     init,
//     //                                 double tol1, double tol2,
//     //                                 size_t max_iter) const;
//     // SurfPoint newton_step(const Point& p, const Vec& w, const Vec& fu,
//     //                       const Vec& fv, const SurfPoint& q) const;
//     // SurfPoint bound_check(const SurfPoint& q) const;

// private:
//     void init_cpoints(const vector<vector<Point>>& p,
//                       const vector<vector<double>>& w);

//     array<size_t, 2> order_;
//     array<vector<double>, 2> knots_;
//     vector<vector<CPoint>> cpoints_;
//     vector<CoxDeBoor<CPoint>> cdb_;
// };

// BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<double>& ku,
//                            const vector<double>& kv,
//                            const vector<vector<Point>>& p,
//                            const vector<vector<double>>& w)
//     : order_{du + 1, dv + 1}
//     , knots_{ku, kv}
//     , cpoints_()
//     , cdb_()
// {
//     init_cpoints(p, w);
// }

// BSplineSurface::Impl::Impl(size_t du, size_t dv, const vector<size_t>&
// ku_mult,
//                            const vector<double>& ku_vals,
//                            const vector<size_t>& kv_mult,
//                            const vector<double>& kv_vals,
//                            const vector<vector<Point>>& p,
//                            const vector<vector<double>>& w)
//     : order_{du + 1, dv + 1}
//     , knots_()
//     , cpoints_()
//     , cdb_()
// {
//     for (size_t i = 0; i < ku_mult.size(); ++i)
//         for (auto j = ku_mult[i]; j > 0; --j)
//             knots_[0].emplace_back(ku_vals[i]);
//     for (size_t i = 0; i < kv_mult.size(); ++i)
//         for (auto j = kv_mult[i]; j > 0; --j)
//             knots_[1].emplace_back(kv_vals[i]);
//     for (auto& i : knots_)
//         i.shrink_to_fit();

//     init_cpoints(p, w);
// }

// Point BSplineSurface::Impl::f(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = cdb_[i].proxy(p.v).get(0);
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Point(cdb_u.proxy(p.u).get(0).p());
// }

// Vec BSplineSurface::Impl::dfu(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = cdb_[i].proxy(p.v).get(0);
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Vec(CPoint::d1(cdb_u.proxy(p.u).range(2)).p());
// }

// Vec BSplineSurface::Impl::dfv(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = CPoint::d1(cdb_[i].proxy(p.v).range(2));
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Vec(cdb_u.proxy(p.u).get(0).p());
// }

// Vec BSplineSurface::Impl::dfuu(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = cdb_[i].proxy(p.v).get(0);
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Vec(CPoint::d2(cdb_u.proxy(p.u).range(3)).p());
// }

// Vec BSplineSurface::Impl::dfvv(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = CPoint::d2(cdb_[i].proxy(p.v).range(3));
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Vec(cdb_u.proxy(p.u).get(0).p());
// }

// Vec BSplineSurface::Impl::dfuv(const SurfPoint& p) const
// {
//     auto n = cpoints_.size();
//     vector<CPoint> cp(n);

//     for (size_t i = 0; i < n; ++i) {
//         cp[i] = CPoint::d1(cdb_[i].proxy(p.v).range(2));
//     }

//     auto cdb_u = CoxDeBoor<CPoint>(order_[0], knots_[0], cp);
//     return Vec(CPoint::d1(cdb_u.proxy(p.u).range(2)).p());
// }

// ostream& BSplineSurface::Impl::print(ostream& os) const
// {
//     fmt::print(os,
//                "{{ \"type\": \"bspline\", \"du\": {0}, \"dv\": {1}, \"ku\":
//                "
//                "{2}, \"kv\": {3}, \"cpoints\": {4} }}",
//                order_[0] - 1, order_[1] - 1, knots_[0], knots_[1],
//                cpoints_);
//     return os;
// }

// SurfPoint BSplineSurface::Impl::project(const Point& p,
//                                         const BSplineSurface& c) const
// {
//     return SurfPoint();
// }

// void BSplineSurface::Impl::init_cpoints(const vector<vector<Point>>& p,
//                                         const vector<vector<double>>& w)
// {
//     auto n = p.size();
//     auto m = p[0].size();

//     cdb_.resize(n);
//     cpoints_.resize(n);
//     if (w.empty()) {
//         for (size_t i = 0; i < n; ++i) {
//             cpoints_[i].resize(m);
//             for (size_t j = 0; j < m; ++j) {
//                 cpoints_[i][j] = CPoint(p[i][j].raw(), 1);
//             }
//             cdb_[i] = CoxDeBoor<CPoint>(order_[1], knots_[1], cpoints_[i]);
//         }

//     } else {
//         for (size_t i = 0; i < n; ++i) {
//             cpoints_[i].resize(m);
//             for (size_t j = 0; j < m; ++j) {
//                 cpoints_[i][j] = CPoint(p[i][j].raw(), w[i][j]);
//             }
//             cdb_[i] = CoxDeBoor<CPoint>(order_[1], knots_[1], cpoints_[i]);
//         }
//     }
// }

// #define I_(k) ((k) % 4)

// SurfPoint BSplineSurface::Impl::project(const Point& p,
//                                         const BSplineSurface& c) const
// {
//     static constexpr array<double, 3> tol = {1e-1, 1e-2, 1e-5};
//     static constexpr size_t dcoeff = 5;
//     static constexpr size_t max_iter = 10;

//     optional<SurfPoint> result;
//     auto min_dist = numeric_limits<double>::max();

//     auto ubegin = knots_[0].front();
//     auto vbegin = knots_[1].front();
//     auto uend = knots_[0].back();
//     auto vend = knots_[1].back();

//     array<SurfPoint, 4> s
//         = {SurfPoint(ubegin, vbegin), SurfPoint(uend, vbegin),
//            SurfPoint(uend, vend), SurfPoint(ubegin, vend)};
//     array<SurfPoint, 4> r;
//     array<double, 4> d;
//     size_t i;

//     auto [k, flag] = is_init(p, s);
//     DEBUG_FMT_("s: {}", s);
//     for (i = 0; !flag && i < max_iter && min(abs(s[0] - s[2])) > tol[1];
//     ++i) {
//         double rd = numeric_limits<double>::max();
//         size_t id = 0;
//         for (size_t j = 0; j < 4; ++j) {
//             r[j] = (s[j] + s[I_(j + 1)]) / 2;
//             d[j] = dist(f(r[j]), p);
//         }

//         auto da = max(d[k], d[I_(k + 2)]);
//         auto db = max(d[I_(k + 1)], d[I_(k + 3)]);
//         if (da < db) {
//             s = {s[k], r[k], r[I_(k + 2)], s[I_(k + 3)]};
//         } else {
//             s = {s[k], s[I_(k + 1)], r[I_(k + 1)], r[I_(k + 3)]};
//         }
//         DEBUG_FMT_("s: {}", s);
//         tie(k, flag) = is_init(p, s);
//     }

//     auto sd = abs(s[0] - s[2]);
//     if (sd.u > tol[0] && sd.v < tol[0]) {
//         DEBUG_FMT_("u project u: [{} {}] v: {}", s[0].u, s[2].u, s[k].v);
//         result = project_u(p, s[k].v, s[0].u, s[2].u, tol[2], tol[2],
//         max_iter,
//                            dcoeff, c);
//     } else if (sd.v > tol[0] && sd.u < tol[0]) {
//         DEBUG_FMT_("v project v: [{} {}] u: {}", s[0].v, s[2].v, s[k].u);
//         result = project_v(p, s[k].u, s[0].v, s[2].v, tol[2], tol[2],
//         max_iter,
//                            dcoeff, c);
//     } else {
//         DEBUG_FMT_("init: {} i: {} flag: {}", s[k], i, flag);
//         result = newton_iter(p, s[k], tol[2], tol[2], max_iter);
//     }

//     if (!result.has_value()) {
//         THROW_(runtime_error, "unable to project point {} on surface {}", p,
//                c);
//     }

//     return result.value();
// }

// pair<size_t, bool>
// BSplineSurface::Impl::is_init(const Point& p,
//                               const array<SurfPoint, 4>& control) const
// {
//     array<Point, 4> q;
//     transform(begin(control), end(control), begin(q),
//               [&](auto& i) { return f(i); });

//     size_t imin = 0;
//     auto min_dist = dist(q[0], p);
//     for (size_t i = 1; i < 4; ++i) {
//         if (auto cur_dist = dist(q[i], p); cur_dist < min_dist) {
//             min_dist = cur_dist;
//             imin = i;
//         }
//     }

//     bool flag = true;
//     auto w = q[imin] - p;
//     for (size_t i = 0; i < 4; ++i) {
//         if (i == imin) {
//             continue;
//         }
//         if (auto dotw = dot(q[i] - q[imin], w); dotw <= 0) {
//             flag = false;
//             break;
//         }
//     }

//     return {size_t(imin), flag};
// }

// optional<SurfPoint> BSplineSurface::Impl::project_u(
//     const Point& p, double v, double ubegin, double uend, double tol1,
//     double tol2, size_t max_iter, size_t dcoeff, const BSplineSurface& c)
//     const
// {
//     if (ubegin >= uend) {
//         swap(ubegin, uend);
//     }

//     optional<SurfPoint> result;
//     auto min_dist = numeric_limits<double>::max();
//     auto udelta = (uend - ubegin) / dcoeff;
//     auto q = SurfPoint(ubegin, v);
//     double du = 0;

//     while (!isnear(q.u, uend, tol2)) {
//         du = min(uend - q.u, udelta * ustep(q, c));
//         DEBUG_FMT_("init: {}", q);
//         auto op = newton_iter(p, q, tol1, tol2, max_iter);
//         if (op.has_value()) {
//             auto s = op.value();
//             if (auto cur_dist = dist(f(s), p); cur_dist < min_dist) {
//                 min_dist = cur_dist;
//                 result = s;
//                 DEBUG_FMT_("cur_dist: {} result {}", cur_dist, s);
//             }
//         }
//         q.u += du;
//     }

//     return result;
// }

// optional<SurfPoint> BSplineSurface::Impl::project_v(
//     const Point& p, double u, double vbegin, double vend, double tol1,
//     double tol2, size_t max_iter, size_t dcoeff, const BSplineSurface& c)
//     const
// {
//     if (vbegin >= vend) {
//         swap(vbegin, vend);
//     }

//     optional<SurfPoint> result;
//     auto min_dist = numeric_limits<double>::max();
//     auto vdelta = (vend - vbegin) / dcoeff;
//     auto q = SurfPoint(u, vbegin);
//     double dv = 0;

//     while (!isnear(q.v, vend, tol2)) {
//         dv = min(vend - q.v, vdelta * vstep(q, c));
//         DEBUG_FMT_("init: {}", q);
//         auto op = newton_iter(p, q, tol1, tol2, max_iter);
//         if (op.has_value()) {
//             auto s = op.value();
//             if (auto cur_dist = dist(f(s), p); cur_dist < min_dist) {
//                 min_dist = cur_dist;
//                 result = s;
//                 DEBUG_FMT_("cur_dist: {} result {}", cur_dist, s);
//             }
//         }
//         q.v += dv;
//     }

//     return result;
// }

// SurfPoint BSplineSurface::Impl::project_s(const Point& p,
//                                           const BSplineSurface& c) const
// {
//     static constexpr auto tol1 = 1e-5;
//     static constexpr auto tol2 = 1e-7;
//     static constexpr auto stol = 1e-10;
//     static constexpr size_t dcoeff = 10;
//     static constexpr size_t max_iter = 100;

//     optional<SurfPoint> result;
//     auto min_dist = numeric_limits<double>::max();

//     double dv = 0;
//     double du = 0;

//     auto ubegin = knots_[0].front();
//     auto vbegin = knots_[1].front();
//     auto uend = knots_[0].back();
//     auto vend = knots_[1].back();

//     auto udelta = (uend - ubegin) / dcoeff;
//     auto vdelta = (vend - vbegin) / dcoeff;

//     SurfPoint q(ubegin, vbegin);
//     Vec n, duu, dvv;

//     while (!isnear(q.u, uend, stol)) {
//         du = min(uend - q.u, udelta * ustep(q, c));

//         q.v = vbegin;
//         while (!isnear(q.v, vend, stol)) {
//             dv = min(vend - q.v, vdelta * vstep(q, c));

//             auto r = SurfPoint({q.u + du, q.v + dv});
//             if (init_between(p, q, r)) {
//                 for (auto& t : {q, r, (q + r) / 2}) {
//                     DEBUG_FMT_("init: {}", t);
//                     auto op = newton_iter(p, t, tol1, tol2, max_iter);
//                     if (op.has_value()) {
//                         auto s = op.value();
//                         if (auto cur_dist = dist(f(s), p);
//                             cur_dist < min_dist) {
//                             min_dist = cur_dist;
//                             result = s;
//                             DEBUG_FMT_("cur_dist: {} result {}", cur_dist,
//                             s);
//                         }
//                     }
//                 }
//             }
//             q.v += dv;
//         }
//         q.u += du;
//     }

//     if (!result.has_value()) {
//         THROW_(runtime_error, "unable to project point {} on surface {}", p,
//                c);
//     }

//     return result.value();
// }

// inline double BSplineSurface::Impl::ustep(const SurfPoint& p,
//                                           const BSplineSurface& c) const
// {
//     double cmul = 1;
//     if (order_[0] > 2) {
//         cmul = norm(dfu(p)) / fabs(dot(c.unit_normal(p), dfuu(p)));
//     }
//     return cmul;
// }

// inline double BSplineSurface::Impl::vstep(const SurfPoint& p,
//                                           const BSplineSurface& c) const
// {
//     double cmul = 1;
//     if (order_[1] > 2) {
//         cmul = norm(dfv(p)) / fabs(dot(c.unit_normal(p), dfvv(p)));
//     }
//     return cmul;
// }

// inline double BSplineSurface::Impl::init_between(const Point& p,
//                                                  const SurfPoint& a,
//                                                  const SurfPoint& b) const
// {
//     auto c = (a + b) / 2;
//     auto d = abs(b - a) / 2;

//     auto fu = dfu(c), fv = dfv(c);
//     auto a1 = dot(fu, f({c.u - d.u, c.v}) - p),
//          a2 = dot(fu, f({c.u + d.u, c.v}) - p),
//          b1 = dot(fv, f({c.u, c.v - d.v}) - p),
//          b2 = dot(fv, f({c.u, c.v + d.v}) - p);

//     return (a1 * a2 <= 0) && (b1 * b2 <= 0);
// }

// optional<SurfPoint> BSplineSurface::Impl::newton_iter(const Point& p,
//                                                       const SurfPoint& init,
//                                                       double tol1, double
//                                                       tol2, size_t max_iter)
//                                                       const
// {
//     optional<SurfPoint> result = nullopt;
//     size_t i;
//     auto u = init;

//     for (i = 0; i < max_iter; ++i) {
//         auto w = Vec(p, f(u));
//         auto fu = dfu(u);
//         auto fv = dfv(u);

//         auto xaa = norm(w);
//         auto xbb = fabs(1 - cos(w, cross(fu, fv)));

//         auto xa = xaa < tol1;
//         auto xb = xbb < tol2;

//         DEBUG_FMT_("norm(w): {} |cos(m, w)|: {}", xaa, xbb);

//         if (xa || xb) {
//             DEBUG_FMT_("xa: {} xb: {}", xa, xb);
//             result = u;
//             break;
//         }

//         auto v = bound_check(u + newton_step(p, w, fu, fv, u));
//         auto xc = norm(fu * (v.u - u.u) + fv * (v.v - u.v)) < tol1;
//         if (xc) {
//             DEBUG_FMT_("xc: {}", xc);
//             result = u;
//             break;
//         }

//         u = v;
//     }

//     if (result.has_value()) {
//         DEBUG_FMT_("converged {}/{}", i, max_iter);
//     } else {
//         DEBUG_FMT_("not converged");
//     }
//     return result;
// }

// inline SurfPoint BSplineSurface::Impl::newton_step(const Point& p,
//                                                    const Vec& w, const Vec&
//                                                    fu, const Vec& fv, const
//                                                    SurfPoint& q) const
// {

//     auto a11 = dot(dfuu(q), w) + sqr(fu);
//     auto a12 = dot(dfuv(q), w) + dot(fu, fv);
//     auto a22 = dot(dfvv(q), w) + sqr(fv);
//     auto b1 = -dot(fu, w);
//     auto b2 = -dot(fv, w);

//     auto d0 = a11 * a22 - sqr(a12);
//     auto d1 = b1 * a22 - b2 * a12;
//     auto d2 = b2 * a11 - b1 * a12;

//     return {d1 / d0, d2 / d0};
// }

// #define R_(i) (knots_[(i)].back())
// #define L_(i) (knots_[(i)].front())
// #define MR_(i, x) (2 * R_(i) - (x))
// #define ML_(i, x) (2 * L_(i) - (x))

// inline SurfPoint BSplineSurface::Impl::bound_check(const SurfPoint& q) const
// {
//     auto r = q;

//     while (r.u < L_(0) || r.u > R_(0)) {
//         if (r.u < L_(0)) {
//             r.u = ML_(0, r.u);
//         }
//         if (r.u > R_(0)) {
//             r.u = MR_(0, r.u);
//         }
//     }

//     while (r.v < L_(1) || r.v > R_(1)) {
//         if (r.v < L_(1)) {
//             r.v = ML_(1, r.v);
//         }
//         if (r.v > R_(1)) {
//             r.v = MR_(1, r.v);
//         }
//     }

//     DEBUG_FMT_("bound check q: {} r: {}", q, r);

//     return r;
// }

// #undef ML_
// #undef MR_
// #undef L_
// #undef R_
