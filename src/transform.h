#pragma once

#include "matrix.h"
#include "vector.h"
#include "point.h"
#include "ray.h"
#include "intersection.h"

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
};

Transform inverse(const Transform &t);
Transform translate(const Vector3f &delta);
Transform scale(float x, float y, float z);

Transform perspective(float vfov, float aspect, float n, float f);
Transform look_at(const Vector3f &eye, const Vector3f &target, const Vector3f &up);

template<typename T>
inline Vector3<T> Transform::operator*(const Vector3<T> &v) const
{
    T x = v.x, y = v.y, z = v.z;
    return Vector3<T>(m.m[0][0] * x + m.m[1][0] * y + m.m[2][0] * z,
                      m.m[0][1] * x + m.m[1][1] * y + m.m[2][1] * z,
                      m.m[0][2] * x + m.m[1][2] * y + m.m[2][2] * z);
}

template<typename T>
inline Point3<T> Transform::operator*(const Point3<T> &p) const
{
    // NOTE: assuming p.w = 1.
    T x = p.x, y = p.y, z = p.z;
    T xp = m.m[0][0] * x + m.m[1][0] * y + m.m[2][0] * z + m.m[3][0];
    T yp = m.m[0][1] * x + m.m[1][1] * y + m.m[2][1] * z + m.m[3][1];
    T zp = m.m[0][2] * x + m.m[1][2] * y + m.m[2][2] * z + m.m[3][2];
    T wp = m.m[0][3] * x + m.m[1][3] * y + m.m[2][3] * z + m.m[3][3];
    if (wp == 1)
        return Point3<T>(xp, yp, zp);
    else
        return Point3<T>(xp, yp, zp) / wp;
}

template<typename T>
inline Normal3<T> Transform::operator*(const Normal3<T> &n) const
{
    T x = n.x, y = n.y, z = n.z;
    return Normal3<T>(inv.m[0][0] * x + inv.m[0][1] * y + inv.m[0][2] * z,
                      inv.m[1][0] * x + inv.m[1][1] * y + inv.m[1][2] * z,
                      inv.m[2][0] * x + inv.m[2][1] * y + inv.m[2][2] * z);
}
