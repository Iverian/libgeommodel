#ifndef GEOM_MODEL_INCLUDE_DOT_H
#define GEOM_MODEL_INCLUDE_DOT_H

template <class U, class V>
double dot(const U& u, const V& v)
{
    double result = 0.;
    auto i = u.cbegin();
    auto j = v.cbegin();
    for (; i != u.cend() && j != v.cend(); ++i, ++j)
        result += (*i) * (*j);
    return result;
};

#endif // GEOM_MODEL_INCLUDE_DOT_H
