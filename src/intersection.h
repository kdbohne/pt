#pragma once

#include "point.h"
#include "normal.h"
#include "ray.h"

struct Entity;

struct Intersection
{
    Point3f p;
    Normal3f n;

    float time;

    const Entity *entity;

    Intersection() {}
    Intersection(const Point3f &p);

    Ray spawn_ray_to(const Intersection &p1) const;
};
