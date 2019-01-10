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
static constexpr float INV_2_PI = 0.15915494309189533577;
static constexpr float PI_OVER_2 = 1.57079632679489661923;
static constexpr float PI_OVER_4 = 0.78539816339744830961;

// TODO: move to random.h when it is created?
static const float ONE_MINUS_EPSILON = 0x1.fffffep-1;

inline int log2int(uint32_t v)
{
    return 31 - __builtin_clz(v);
}

template<typename T>
inline bool is_power_of_2(T v)
{
    return v && !(v & (v - 1));
}

inline int32_t round_up_power_of_2(int32_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return v + 1;
}

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

inline float lerp(float t, float v0, float v1)
{
    return (1 - t) * v0 + t * v1;
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

template<typename Predicate>
inline int find_interval(int size, const Predicate &pred)
{
    int first = 0;
    int len = size;
    while (len > 0)
    {
        int half = len >> 1;
        int middle = first + half;
        if (pred(middle))
        {
            first = middle + 1;
            len -= half + 1;
        }
        else
        {
            len = half;
        }
    }

    return clamp(first - 1, 0, size - 2);
}
