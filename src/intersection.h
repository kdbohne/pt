#pragma once

#include "point.h"
#include "vector.h"
#include "ray.h"

struct Entity;

struct Intersection
{
    Point3f p;
    Vector3f n;
    Vector3f wo;

    float time;

    const Entity *entity;

    Intersection() {}
    Intersection(const Point3f &p, const Vector3f &n, const Vector3f &wo, float time)
        : p(p), n(n), wo(wo), time(time), entity(nullptr)
    {
    }

    Ray spawn_ray(const Vector3f &d) const
    {
        return Ray(p, d, INFINITY, time);
    }
};
