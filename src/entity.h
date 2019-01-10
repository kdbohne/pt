#pragma once

#include "spectrum.h"
#include "point.h"
#include "vector.h"

struct Geometry;
struct AreaLight;
struct Material;
struct Ray;
struct Intersection;

struct Entity
{
    Geometry *geometry;
    AreaLight *area_light;
    Material *material;

    Entity() {}
    Entity(Geometry *geometry, AreaLight *light, Material *material)
        : geometry(geometry), area_light(light), material(material) {}

    bool intersect(const Ray &ray, Intersection *intersection) const;
};
