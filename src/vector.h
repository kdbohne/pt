#pragma once

#include <cmath>
#include "common.h"

template<typename T>
struct Normal3;

template<typename T>
struct Vector2
{
    T x, y;

    Vector2() : x(0), y(0) {}
    Vector2(T x, T y) : x(x), y(y) {}

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

    Vector2<T> operator/(T s) const
    {
        return Vector2<T>(x / s, y / s);
    }
};

typedef Vector2<float> Vector2f;
typedef Vector2<int>   Vector2i;

template<typename T>
struct Vector3
{
    T x, y, z;

    Vector3() : x(0), y(0), z(0) {}
    Vector3(T s) : x(s), y(s), z(s) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    explicit Vector3(const Normal3<T> &n) : x(n.x), y(n.y), z(n.z) {}

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

    Vector3<T> operator+(const Vector3<T> &v) const
    {
        return Vector3<T>(x + v.x, y + v.y, z + v.z);
    }

    Vector3<T> operator-(const Vector3<T> &v) const
    {
        return Vector3<T>(x - v.x, y - v.y, z - v.z);
    }

    Vector3<T> &operator+=(const Vector3<T> &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vector3<T> &operator-=(const Vector3<T> &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    Vector3<T> operator*(T s) const
    {
        return Vector3<T>(x * s, y * s, z * s);
    }

    Vector3<T> operator/(T s) const
    {
        return Vector3<T>(x / s, y / s, z / s);
    }

    Vector3<T> operator-() const
    {
        return Vector3<T>(-x, -y, -z);
    }

    float length_squared() const
    {
        return (x * x) + (y * y) + (z * z);
    }

    float length() const
    {
        return std::sqrt(length_squared());
    }
};

typedef Vector3<float> Vector3f;
typedef Vector3<int>   Vector3i;

template<typename T>
inline Vector3<T> operator*(T s, const Vector3<T> &v)
{
    return v * s;
}

template<typename T>
inline Vector3<T> normalize(const Vector3<T> &v)
{
    return v / v.length();
}

template<typename T>
inline Vector3<T> abs(const Vector3<T> &v)
{
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template<typename T>
inline T dot(const Vector3<T> &a, const Vector3<T> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template<typename T>
inline T dot(const Vector3<T> &a, const Normal3<T> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template<typename T>
inline T dot(const Normal3<T> &a, const Vector3<T> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template<typename T>
inline T abs_dot(const Vector3<T> &a, const Vector3<T> &b)
{
    return std::abs(dot(a, b));
}

template<typename T>
inline T abs_dot(const Vector3<T> &a, const Normal3<T> &b)
{
    return std::abs(dot(a, b));
}

template<typename T>
inline T abs_dot(const Normal3<T> &a, const Vector3<T> &b)
{
    return std::abs(dot(a, b));
}

template<typename T>
inline Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b)
{
    double ax = a.x, ay = a.y, az = a.z;
    double bx = b.x, by = b.y, bz = b.z;
    return Vector3<T>((ay * bz) - (az * by),
                      (az * bx) - (ax * bz),
                      (ax * by) - (ay * bx));
}

template<typename T>
inline Vector3<T> cross(const Normal3<T> &a, const Vector3<T> &b)
{
    double ax = a.x, ay = a.y, az = a.z;
    double bx = b.x, by = b.y, bz = b.z;
    return Vector3<T>((ay * bz) - (az * by),
                      (az * bx) - (ax * bz),
                      (ax * by) - (ay * bx));
}

template<typename T>
inline int max_dimension(const Vector3<T> &v)
{
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template<typename T>
inline Vector3<T> permute(const Vector3<T> &v, int x, int y, int z)
{
    return Vector3<T>(v[x], v[y], v[z]);
}

template<typename T>
inline void coordinate_system(const Vector3<T> &v0, Vector3<T> *v1, Vector3<T> *v2)
{
    // TODO: does this work for right-handed coordinate systems?
    if (std::abs(v0.x) > std::abs(v0.y))
        *v1 = Vector3<T>(-v0.z, 0, v0.x) / std::sqrt(v0.x * v0.x + v0.z * v0.z);
    else
        *v1 = Vector3<T>(0, v0.z, -v0.y) / std::sqrt(v0.y * v0.y + v0.z * v0.z);

    *v2 = cross(v0, *v1);
}

inline Vector3f spherical_direction(float sin_theta, float cos_theta, float phi)
{
    return Vector3f(sin_theta * std::cos(phi),
                    sin_theta * std::sin(phi),
                    cos_theta);
}

inline Vector3f spherical_direction(float sin_theta, float cos_theta, float phi,
                                    const Vector3f &x, const Vector3f &y, const Vector3f &z)
{
    return sin_theta * std::cos(phi) * x +
           sin_theta * std::sin(phi) * y +
           cos_theta * z;
}
