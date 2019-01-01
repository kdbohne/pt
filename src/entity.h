#pragma once

#include "spectrum.h"
#include "point.h"
#include "vector.h"

struct Geometry;
struct AreaLight;
struct Bsdf;
struct Ray;
struct Intersection;

struct Entity
{
    Geometry *geometry;
    AreaLight *area_light;
    Bsdf *bsdf;

    Entity() {}
    Entity(Geometry *geometry, AreaLight *light, Bsdf *bsdf)
        : geometry(geometry), area_light(light), bsdf(bsdf) {}

    bool intersect(const Ray &ray, Intersection *intersection) const;
};
