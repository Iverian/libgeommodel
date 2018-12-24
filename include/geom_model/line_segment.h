#ifndef GEOM_MODEL_INCLUDE_LINE_SEGMENT_H_
#define GEOM_MODEL_INCLUDE_LINE_SEGMENT_H_

#include "oriented_edge.h"

class LineSegment : OrientedEdge {
public:
    LineSegment(const Point& a, const Point& b);

private:
};

#endif // GEOM_MODEL_INCLUDE_LINE_SEGMENT_H_