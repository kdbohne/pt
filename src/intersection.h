#pragma once

#include "point.h"
#include "normal.h"
#include "ray.h"

struct Entity;

struct Intersection
{
    Point3f p;
    Normal3f n;
    Vector3f wo;

    Vector3f t, b;

    float time;

    const Entity *entity;

    Intersection() {}
    Intersection(const Point3f &p);

    Ray spawn_ray(const Vector3f &d) const;
    Ray spawn_ray_to(const Intersection &p1) const;
};
