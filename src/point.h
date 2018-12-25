#pragma once

#include "vector.h"

template<typename T>
struct Point2
{
    T x, y;

    Point2() : x(0), y(0) {}
    Point2(T x, T y) : x(x), y(y) {}

    template<typename U>
    explicit Point2(const Vector2<U> &v) : x((T)v.x), y((T)v.y) {}

    template<typename U>
    explicit operator Vector2<U>() const
    {
        return Vector2<U>((U)x, (U)y);
    }

    T operator[](int i) const
    {
        assert((i >= 0) && (i <= 1));
        return *(&x + i);
    }

    T &operator[](int i)
    {
        assert((i >= 0) && (i <= 1));
        return *(&x + i);
    }

    Point2<T> operator+(const Point2<T> &p) const
    {
        return Point2<T>(x + p.x, y + p.y);
    }

    Point2<T> operator+(const Vector2<T> &v) const
    {
        return Point2<T>(x + v.x, y + v.y);
    }

    Vector2<T> operator-(const Point2<T> &p) const
    {
        return Vector2<T>(x - p.x, y - p.y);
    }

    Point2<T> operator-(const Vector2<T> &v) const
    {
        return Point2<T>(x - v.x, y - v.y);
    }

    Point2<T> operator*(T s) const
    {
        return Point2<T>(x * s, y * s);
    }

    Point2<T> operator/(T s) const
    {
        return Point2<T>(x / s, y / s);
    }
};

typedef Point2<float> Point2f;
typedef Point2<int>   Point2i;

template<typename T>
inline Point2<T> operator*(T s, const Point2<T> &p)
{
    return p * s;
}

template<typename T>
struct Point3
{
    T x, y, z;

    Point3() : x(0), y(0), z(0) {}
    Point3(T x, T y, T z) : x(x), y(y), z(z) {}

    template<typename U>
    explicit Point3(const Vector3<U> &v) : x((T)v.x), y((T)v.y), z((T)v.z) {}

    template<typename U>
    explicit operator Vector3<U>() const
    {
        return Vector3<U>((U)x, (U)y, (U)z);
    }

    T operator[](int i) const
    {
        assert((i >= 0) && (i <= 2));
        return *(&x + i);
    }

    T &operator[](int i)
    {
        assert((i >= 0) && (i <= 2));
        return *(&x + i);
    }

    Point3<T> operator+(const Point3<T> &p) const
    {
        return Point3<T>(x + p.x, y + p.y, z + p.z);
    }

    Point3<T> operator+(const Vector3<T> &v) const
    {
        return Point3<T>(x + v.x, y + v.y, z + v.z);
    }

    Vector3<T> operator-(const Point3<T> &p) const
    {
        return Vector3<T>(x - p.x, y - p.y, z - p.z);
    }

    Point3<T> operator*(T s) const
    {
        return Point3<T>(x * s, y * s, z * s);
    }

    Point3<T> &operator*=(T s)
    {
        x *= s;
        y *= s;
        z *= s;
        return *this;
    }

    Point3<T> operator/(T s) const
    {
        return Point3<T>(x / s, y / s, z / s);
    }
};

typedef Point3<float> Point3f;
typedef Point3<int>   Point3i;

template<typename T>
inline Point3<T> operator*(T s, const Point3<T> &p)
{
    return p * s;
}

template<typename T>
inline float distance(const Point3<T> &a, const Point3<T> &b)
{
    return (a - b).length();
}

template<typename T>
inline float distance_squared(const Point3<T> &a, const Point3<T> &b)
{
    return (a - b).length_squared();
}
