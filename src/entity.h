#pragma once

#include "spectrum.h"
#include "point.h"
#include "vector.h"

struct Geometry;
struct Light;
struct Bsdf;
struct Ray;
struct Intersection;

struct Entity
{
    Geometry *geometry;
    Light *light;
    Bsdf *bsdf;

    Entity() {}
    Entity(Geometry *g, Light *l, Bsdf *b) : geometry(g), light(l), bsdf(b) {}

    bool intersect(const Ray &ray, Intersection *intersection) const;
};
