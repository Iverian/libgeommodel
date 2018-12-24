#include <geom_model/line.h>
#include <geom_model/line_segment.h>

using namespace std;

LineSegment::LineSegment(const Point& a, const Point& b)
    : OrientedEdge(Edge(make_shared<Line>(Vec(b - a), a), 0, 1), true)
{
}