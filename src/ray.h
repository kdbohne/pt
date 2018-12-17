#pragma once

#include "vector.h"
#include "point.h"

struct Ray
{
    Point3f o;
    Vector3f d;
    mutable float tmax;

    float time;

    Ray() : o(Point3f()), d(Vector3f()), tmax(INFINITY), time(0) {}
    Ray(const Point3f &o, const Vector3f &d, float tmax, float time) : o(o), d(d), tmax(tmax), time(time) {}

    Point3f evaluate(float t) const
    {
        return o + t * d;
    }
};
