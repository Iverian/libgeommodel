#ifndef GEOM_MODEL_INCLUDE_EDGE_H_
#define GEOM_MODEL_INCLUDE_EDGE_H_

#include "abstract_curve.h"

#include <functional>
#include <iostream>
#include <memory>

class Edge;

std::ostream& operator<<(std::ostream& os, const Edge& edge);

class Edge {
    struct Impl;

public:
    Edge();
    Edge(const std::shared_ptr<AbstractCurve>& curve, const Point& begin,
         const Point& end);
    Edge(const std::shared_ptr<AbstractCurve>& curve, double begin,
         double end);

    Point f(double t) const;
    Point front() const;
    Point back() const;
    double pfront() const;
    double pback() const;
    double length() const;

    const AbstractCurve& curve() const;

    friend bool operator==(const Edge& lhs, const Edge& rhs);
    friend bool operator!=(const Edge& lhs, const Edge& rhs);

private:
    std::shared_ptr<Impl> pimpl_;
};

namespace std {
template <>
struct hash<Edge> {
    size_t operator()(const Edge& key) const
    {
        return (chasher_(key.curve().shared_from_this()) << 2)
            ^ (dhasher_(key.pfront()) << 1) ^ (dhasher_(key.pback()));
    }

private:
    hash<shared_ptr<const AbstractCurve>> chasher_;
    hash<double> dhasher_;
};
}

#endif // GEOM_MODEL_INCLUDE_EDGE_H_
