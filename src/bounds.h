#pragma once

#include "point.h"

template<typename T>
struct Bounds2
{
    Point2<T> min, max;

    Bounds2()
    {
        T min_value = std::numeric_limits<T>::lowest();
        T max_value = std::numeric_limits<T>::max();

        min = Point2<T>(max_value, max_value);
        max = Point2<T>(min_value, min_value);
    }

    Bounds2(const Point2<T> &p0, const Point2<T> &p1)
        : min(std::min(p0.x, p1.x), std::min(p0.y, p1.y)),
          max(std::max(p0.x, p1.x), std::max(p0.y, p1.y))
    {
    }
};

typedef Bounds2<float> Bounds2f;

template<typename T>
struct Bounds3
{
    Point3<T> min, max;

    Bounds3()
    {
        T min_value = std::numeric_limits<T>::lowest();
        T max_value = std::numeric_limits<T>::max();

        min = Point3<T>(max_value, max_value, max_value);
        max = Point3<T>(min_value, min_value, min_value);
    }

    Bounds3(const Point3<T> &p0, const Point3<T> &p1)
        : min(std::min(p0.x, p1.x), std::min(p0.y, p1.y), std::min(p0.z, p1.z)),
          max(std::max(p0.x, p1.x), std::max(p0.y, p1.y), std::max(p0.z, p1.z))
    {
    }

    void bounding_sphere(Point3<T> *center, float *radius)
    {
        *center = (min + max) / 2;
        *radius = inside(*center, *this) ? distance(*center, max) : 0;
    }
};

typedef Bounds3<float> Bounds3f;

template<typename T>
inline Point3<T> min(const Point3<T> &a, const Point3<T> &b)
{
    return Point3<T>(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

template<typename T>
inline Point3<T> max(const Point3<T> &a, const Point3<T> &b)
{
    return Point3<T>(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

// TODO: different name? have to capitalize "union" to avoid clashing with keyword
template<typename T>
inline Bounds3<T> Union(const Bounds3<T> &a, const Bounds3<T> &b)
{
    Bounds3<T> result;
    result.min = min(a.min, b.min);
    result.max = max(a.max, b.max);

    return result;
}

// TODO: different name? have to capitalize "union" to avoid clashing with keyword
template<typename T>
inline Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p)
{
    Bounds3<T> result;
    result.min = min(b.min, p);
    result.max = max(b.max, p);

    return result;
}

template<typename T>
inline bool inside(const Point3<T> &p, const Bounds3<T> &b)
{
    return (((p.x >= b.min.x) && (p.y >= b.min.y) && (p.z >= b.min.z)) &&
            ((p.x <= b.max.x) && (p.y <= b.max.y) && (p.z <= b.max.z)));
}
