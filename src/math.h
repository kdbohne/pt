#pragma once

#include <limits>
#include <cmath>

// TODO: is this needed?
#undef min
#undef max

#undef INFINITY
static constexpr float INFINITY = std::numeric_limits<float>::infinity();
static constexpr float PI = 3.14159265358979323846;
static constexpr float INV_PI = 0.31830988618379067154;
static constexpr float PI_OVER_2 = 1.57079632679489661923;
static constexpr float PI_OVER_4 = 0.78539816339744830961;

inline float radians(float deg)
{
    return deg * (PI / 180.0f);
}

template<typename T, typename U, typename V>
inline T clamp(T v, U min, V max)
{
    if (v < min)
        return min;
    if (v > max)
        return max;
    return v;
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
