#include "bspline_curve_impl.h"
#include <util/debug.h>

#include <vector>

namespace gm {

struct DistanceCurve {
public:
    DistanceCurve();
    DistanceCurve(const BezierPatch& patch, const Point& p);

    DistanceCurve& reconfig(DistanceCurve& out, double a) const __GM_NOEXCEPT_RELEASE__;

    size_t order() const noexcept;
    std::vector<double> knots() const noexcept;
    std::vector<DPoint> points() const noexcept;
    double pfront() const noexcept;
    double pback() const noexcept;

    double f(double u) const noexcept;
    double df(double u) const noexcept;
    double df2(double u) const noexcept;

    bool min_between(double a, double b) const noexcept;
    std::optional<double> minimize(double u0) const noexcept;
    std::optional<size_t> peak_point() const noexcept;


protected:
    double bound_check(double u) const noexcept;
    double armijo_step(double u, double h) const noexcept;

private:
    size_t order_;
    std::vector<double> knots_;
    std::vector<DPoint> points_;
    CoxDeBoor<DPoint> cdb_;
};

struct CurveProjector {
public:
    explicit CurveProjector(const BSplineCurve::Impl& impl);
    double get(const Point& p) const;

    bool is_candidate(double a, const DistanceCurve& c) const noexcept;

private:
    std::vector<BezierPatch> patches_;
    Point front_;
    Point back_;
};

} // namespace gm