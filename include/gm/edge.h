#ifndef GEOM_MODEL_INCLUDE_GM_EDGE_H_
#define GEOM_MODEL_INCLUDE_GM_EDGE_H_

#include "abstract_curve.h"

#include <functional>
#include <iostream>
#include <memory>

namespace gm {

class Edge : public AbstractCurve {
public:
    Edge(std::shared_ptr<AbstractCurve> curve, const Point& front,
         const Point& back);
    Edge(std::shared_ptr<AbstractCurve> curve, double pfront, double pback);

    std::shared_ptr<AbstractCurve> curve() const noexcept;
    Point front() const noexcept;
    Point back() const noexcept;
    double pfront() const noexcept;
    double pback() const noexcept;

    Point f(double u) const noexcept override;
    Vec df(double u) const noexcept override;
    Vec df2(double u) const noexcept override;
    double project(const Point& p) const override;
    std::optional<double> project_greater(const Point& p, double min) const
        noexcept override;
    double approx_length(double begin, double end, size_t n) const override;

    friend bool operator==(const Edge& lhs, const Edge& rhs) noexcept;
    friend bool operator!=(const Edge& lhs, const Edge& rhs) noexcept;

protected:
    std::ostream& print(std::ostream& os) const override;

private:
    double param(double t) const noexcept;
    double param_rev(double u) const noexcept;

    double pfront_;
    double pback_;
    std::shared_ptr<AbstractCurve> curve_;
};

} // namespace gm

namespace std {
template <>
struct hash<gm::Edge> {
    size_t operator()(const gm::Edge& key) const
    {
        return (chasher_(key.curve()) << 2) ^ (dhasher_(key.pfront()) << 1)
            ^ (dhasher_(key.pback()));
    }

private:
    hash<shared_ptr<gm::AbstractCurve>> chasher_;
    hash<double> dhasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_GM_EDGE_H_
