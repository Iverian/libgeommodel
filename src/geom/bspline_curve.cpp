
#include <geom_model/bspline_curve.h>
#include <geom_model/point.h>

#include "cox_de_boor.h"
#include <util/math.h>
#include <util/util.h>

#include <algorithm>
#include <optional>
#include <stdexcept>
#include <vector>

#include <fmt/ostream.h>

using namespace std;

struct BSplineCurve::Impl {
    Impl(size_t degree, vector<double> knots, vector<Point> control_points,
         vector<double> weights = vector<double>());
    Impl(size_t degree, const vector<size_t>& knot_mult,
         const vector<double>& knot_list, vector<Point> control_points,
         vector<double> weights = vector<double>());

    Point f(double u) const;
    Vec df(double u) const;
    Vec df2(double u) const;
    ostream& print(ostream& os) const;
    double project(const Point& p, const BSplineCurve& parent) const;

    const double& param_front() const;
    const double& param_back() const;

    size_t get_order() const;
    const vector<double>& get_knots() const;
    const vector<Point>& get_control_points() const;
    const vector<double>& get_weights() const;

    void init_cdb();

private:
    size_t order_;
    vector<double> knots_;
    vector<Point> control_points_;
    vector<double> weights_;
    CoxDeBoor cdb_;
};

BSplineCurve::~BSplineCurve() = default;
BSplineCurve::BSplineCurve(BSplineCurve&&) noexcept = default;
BSplineCurve& BSplineCurve::operator=(BSplineCurve&&) noexcept = default;
BSplineCurve::BSplineCurve(const BSplineCurve& other)
    : pimpl_(make_unique<Impl>(*other.pimpl_))
{
}
BSplineCurve& BSplineCurve::operator=(const BSplineCurve& other)
{
    *pimpl_ = *other.pimpl_;
    return *this;
}

BSplineCurve::BSplineCurve(size_t degree, vector<double> knots,
                           vector<Point> control_points,
                           vector<double> weights)
    : pimpl_(make_unique<BSplineCurve::Impl>(degree, knots, control_points,
                                             weights))
{
}

BSplineCurve::BSplineCurve(size_t degree,
                           const vector<size_t>& knot_multiplies,
                           const vector<double>& knot_list,
                           vector<Point> control_points,
                           vector<double> weights)
    : pimpl_(make_unique<BSplineCurve::Impl>(degree, knot_multiplies,
                                             knot_list, control_points))
{
}

Point BSplineCurve::f(double u) const
{
    return pimpl_->f(u);
}

Vec BSplineCurve::df(double u) const
{
    return pimpl_->df(u);
}

Vec BSplineCurve::df2(double u) const
{
    return pimpl_->df2(u);
}

ostream& BSplineCurve::print(ostream& os) const
{
    return pimpl_->print(os);
}

double BSplineCurve::project(const Point& p) const
{
    return pimpl_->project(p, *this);
}

const double& BSplineCurve::param_front() const
{
    return pimpl_->param_front();
}

const double& BSplineCurve::param_back() const
{
    return pimpl_->param_back();
}

size_t BSplineCurve::get_order() const
{
    return pimpl_->get_order();
}

const vector<double>& BSplineCurve::get_knots() const
{
    return pimpl_->get_knots();
}

const vector<Point>& BSplineCurve::get_control_points() const
{
    return pimpl_->get_control_points();
}

const vector<double>& BSplineCurve::get_weights() const
{
    return pimpl_->get_weights();
}

BSplineCurve::Impl::Impl(size_t degree, vector<double> knots,
                         vector<Point> control_points, vector<double> weights)
    : order_(degree + 1)
    , knots_(move(knots))
    , control_points_(move(control_points))
    , weights_(move(weights))
    , cdb_()
{
    if (weights_.empty()) {
        weights_.resize(control_points_.size(), 1);
    }
    init_cdb();
}

BSplineCurve::Impl::Impl(size_t degree, const vector<size_t>& knot_mult,
                         const vector<double>& knot_list,
                         vector<Point> control_points, vector<double> weights)
    // : cdb_(degree + 1, knot_multiplies, knot_list, control_points, weights)
    : order_(degree + 1)
    , knots_()
    , control_points_(control_points)
    , weights_(weights)
    , cdb_()
{
    if (weights_.empty()) {
        weights_.resize(control_points_.size(), 1);
    }
    for (size_t i = 0; i < knot_mult.size(); ++i) {
        for (size_t j = knot_mult[i]; j > 0; --j) {
            knots_.emplace_back(knot_list[i]);
        }
    }
    init_cdb();
}

Point BSplineCurve::Impl::f(double u) const
{
    return cdb_.get_proxy(u).f();
}

Vec BSplineCurve::Impl::df(double u) const
{
    return cdb_.get_proxy(u).df();
}

Vec BSplineCurve::Impl::df2(double u) const
{
    return cdb_.get_proxy(u).df2();
}

const double& BSplineCurve::Impl::param_front() const
{
    return knots_.front();
}

const double& BSplineCurve::Impl::param_back() const
{
    return knots_.back();
}

ostream& BSplineCurve::Impl::print(ostream& os) const
{
    fmt::print(os,
               "{{ \"type\": \"bspline\", \"deg\": {0}, \"knots\": {1}, "
               "\"points\": {2}, \"weights\": {3} }}",
               order_, knots_, control_points_, weights_);
    return os;
}

#define MESH_ELEM(i) (param_front() + (param_back() * (i)) / mesh_size)

double BSplineCurve::Impl::project(const Point& p,
                                   const BSplineCurve& parent) const
{
    static constexpr size_t mesh_size = 10;

    vector<double> results;
    for (size_t i = 0; i < mesh_size; ++i) {
        auto a = MESH_ELEM(i), b = MESH_ELEM(i + 1);
        if (parent.is_init_in_interval(p, {a, b})) {
            auto projection = parent.project_iterative(p, (a + b) / 2);
            if (projection.has_value()) {
                results.emplace_back(projection.value());
            }
        }
    }
    if (results.empty()) {
        throw runtime_error(
            fmt::format("Unable project point {0} to curve {1}", p, parent));
    }
    return *min_element(cbegin(results), cend(results),
                        [&](const auto& lhs, const auto& rhs) {
                            return dist(f(lhs), p) < dist(f(rhs), p);
                        });
}

#undef MESH_ELEM

size_t BSplineCurve::Impl::get_order() const
{
    return order_;
}

const vector<double>& BSplineCurve::Impl::get_knots() const
{
    return knots_;
}

const vector<Point>& BSplineCurve::Impl::get_control_points() const
{
    return control_points_;
}

const vector<double>& BSplineCurve::Impl::get_weights() const
{
    return weights_;
}

void BSplineCurve::Impl::init_cdb()
{
    cdb_ = CoxDeBoor(order_, knots_, control_points_, weights_);
}
