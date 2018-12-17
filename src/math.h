#pragma once

#include <limits>
#include <cmath>

#undef INFINITY
static constexpr float INFINITY = std::numeric_limits<float>::infinity();
static constexpr float PI = 3.14159265358979323846;

inline float radians(float deg)
{
    return deg * (PI / 180.0f);
}

// Numerically-stable quadratic equation as written in Physically Based
// Rendering, Third Edition, pp. 1079-1080.
inline bool quadratic(float a, float b, float c, float *t0, float *t1)
{
    double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
    if (discrim < 0)
        return false;

    double root_discrim = std::sqrt(discrim);

    double q;
    if (b < 0)
        q = -0.5 * (b - root_discrim);
    else
        q = -0.5 * (b + root_discrim);

    *t0 = q / a;
    *t1 = c / q;

    if (*t0 > *t1)
        std::swap(*t0, *t1);

    return true;
}
