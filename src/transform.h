#pragma once

#include "matrix.h"
#include "vector.h"
#include "point.h"
#include "ray.h"
#include "intersection.h"
#include "bounds.h"

struct Transform
{
    Matrix4x4 m;
    Matrix4x4 inv;

    Transform() {}
    Transform(const Matrix4x4 &m) : m(m), inv(inverse(m)) {}
    Transform(const Matrix4x4 &m, const Matrix4x4 &inv) : m(m), inv(inv) {}

    Transform operator*(const Transform &t) const;
    Ray operator*(const Ray &r) const;
    Intersection operator*(const Intersection &is) const;

    template<typename T>
    inline Vector3<T> operator*(const Vector3<T> &v) const;
    template<typename T>
    inline Point3<T> operator*(const Point3<T> &p) const;
    template<typename T>
    inline Normal3<T> operator*(const Normal3<T> &n) const;
    template<typename T>
    inline Bounds3<T> operator*(const Bounds3<T> &b) const;
};

Transform inverse(const Transform &t);
Transform translate(const Vector3f &delta);
Transform rotate(float angle, const Vector3f &axis);
Transform scale(float x, float y, float z);

Transform perspective(float vfov, float aspect, float n, float f);
Transform look_at(const Point3f &eye, const Point3f &target, const Vector3f &up);

template<typename T>
inline Vector3<T> Transform::operator*(const Vector3<T> &v) const
{
    return Vector3<T>(m.m[0][0] * v.x + m.m[1][0] * v.y + m.m[2][0] * v.z,
                      m.m[0][1] * v.x + m.m[1][1] * v.y + m.m[2][1] * v.z,
                      m.m[0][2] * v.x + m.m[1][2] * v.y + m.m[2][2] * v.z);
}

template<typename T>
inline Point3<T> Transform::operator*(const Point3<T> &p) const
{
    // NOTE: assuming p.w = 1.
    T xp = m.m[0][0] * p.x + m.m[1][0] * p.y + m.m[2][0] * p.z + m.m[3][0];
    T yp = m.m[0][1] * p.x + m.m[1][1] * p.y + m.m[2][1] * p.z + m.m[3][1];
    T zp = m.m[0][2] * p.x + m.m[1][2] * p.y + m.m[2][2] * p.z + m.m[3][2];
    T wp = m.m[0][3] * p.x + m.m[1][3] * p.y + m.m[2][3] * p.z + m.m[3][3];
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template<typename T>
inline Normal3<T> Transform::operator*(const Normal3<T> &n) const
{
    return Normal3<T>(inv.m[0][0] * n.x + inv.m[0][1] * n.y + inv.m[0][2] * n.z,
                      inv.m[1][0] * n.x + inv.m[1][1] * n.y + inv.m[1][2] * n.z,
                      inv.m[2][0] * n.x + inv.m[2][1] * n.y + inv.m[2][2] * n.z);
}

template<typename T>
inline Bounds3<T> Transform::operator*(const Bounds3<T> &b) const
{
    // TODO: optimize?
    Bounds3<T> result;
    result = Union(result, (*this) * Point3<T>(b.min.x, b.min.y, b.min.z));
    result = Union(result, (*this) * Point3<T>(b.min.x, b.min.y, b.max.z));
    result = Union(result, (*this) * Point3<T>(b.min.x, b.max.y, b.min.z));
    result = Union(result, (*this) * Point3<T>(b.min.x, b.max.y, b.max.z));
    result = Union(result, (*this) * Point3<T>(b.max.x, b.min.y, b.min.z));
    result = Union(result, (*this) * Point3<T>(b.max.x, b.min.y, b.max.z));
    result = Union(result, (*this) * Point3<T>(b.max.x, b.max.y, b.min.z));
    result = Union(result, (*this) * Point3<T>(b.max.x, b.max.y, b.max.z));

    return result;
}
