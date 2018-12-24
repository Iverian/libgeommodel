#ifndef GEOM_MODEL_INCLUDE_ORIENTED_EDGE_H_
#define GEOM_MODEL_INCLUDE_ORIENTED_EDGE_H_

#include "edge.h"

class OrientedEdge;

std::ostream& operator<<(std::ostream& os, const OrientedEdge& oedge);

class OrientedEdge {
public:
    ~OrientedEdge();
    OrientedEdge(OrientedEdge&&) noexcept;
    OrientedEdge& operator=(OrientedEdge&&) noexcept;
    OrientedEdge(const OrientedEdge& other);
    OrientedEdge& operator=(const OrientedEdge& other);

    OrientedEdge(const Edge& edge, bool orientation);
    OrientedEdge(Edge&& edge, bool orienation);

    const Edge& edge() const;
    bool orienation() const;

    Point f(double t) const;
    Point front() const;
    Point back() const;
    double pfront() const;
    double pback() const;

    double length() const;

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

#endif // GEOM_MODEL_INCLUDE_ORIENTED_EDGE_H_
