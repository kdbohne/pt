#pragma once

#include <cmath>

template<typename T>
struct Vector3;
template<typename T>
struct Point3;

template<typename T>
struct Normal3
{
    T x, y, z;

    Normal3() : x(0), y(0), z(0) {}
    Normal3(T s) : x(s), y(s), z(s) {}
    Normal3(T x, T y, T z) : x(x), y(y), z(z) {}

    template<typename U>
    explicit Normal3(const Vector3<U> &v) : x((T)v.x), y((T)v.y), z((T)v.z) {}
    template<typename U>
    explicit Normal3(const Point3<U> &p) : x((T)p.x), y((T)p.y), z((T)p.z) {}

    Normal3<T> operator+(const Normal3<T> &n) const
    {
        return Normal3<T>(x + n.x, y + n.y, z + n.z);
    }

    Normal3<T> operator*(T s) const
    {
        return Normal3<T>(x * s, y * s, z * s);
    }

    Normal3<T> operator/(T s) const
    {
        return Normal3<T>(x / s, y / s, z / s);
    }

    Normal3<T> operator-() const
    {
        return Normal3<T>(-x, -y, -z);
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

typedef Normal3<float> Normal3f;

template<typename T>
inline Normal3<T> operator*(float s, const Normal3<T> &n)
{
    return n * s;
}

template<typename T>
inline Normal3<T> normalize(const Normal3<T> &n)
{
    return n / n.length();
}

template<typename T>
inline Normal3<T> face_forward(const Normal3<T> &n, const Vector3f &v)
{
    return (dot(n, v) < 0) ? -n : n;
}
