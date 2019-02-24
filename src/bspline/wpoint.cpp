#include <bspline/wpoint.h>

namespace gm {

Point pget(const CPoint& p) __GM_NOEXCEPT_RELEASE__
{
    auto& w = p.w();
    check_ifd(!cmp::zero(w), "Zero weight");

    Point result;
    for (CPoint::size_type i = 0; i < p.size(); ++i) {
        result[i] = p[i] / w;
    }

    return result;
}

double pget(const DPoint& p) __GM_NOEXCEPT_RELEASE__
{
    auto& w = p.w();
    check_ifd(!cmp::zero(w), "Zero weight");

    return p[0] / w;
}

} // namespace gm
