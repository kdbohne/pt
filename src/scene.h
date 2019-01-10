#pragma once

#include "entity.h"

#include <vector>

struct Light;
struct Ray;
struct Intersection;

struct Scene
{
    std::vector<Entity> entities;
    std::vector<Light *> lights;

    bool intersect(const Ray &ray, Intersection *intersection) const;
};
