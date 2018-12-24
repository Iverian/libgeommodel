#include <geom_model/parametric_point.h>

#include <fmt/ostream.h>

using namespace std;

ostream& operator<<(ostream& os, const ParametricPoint& p)
{
    fmt::print(os, "[{0}, {1}]", p.u, p.v);
    return os;
}
