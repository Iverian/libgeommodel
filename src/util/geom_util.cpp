#include <cmath>
#include <geom_model/geom_util.h>

using namespace std;

bool isnear(double lhs, double rhs, double eps)
{
    return fabs(rhs - lhs) < eps;
}
