#include "scene.h"
#include "geometry.h"

bool Scene::intersect(const Ray &ray, Intersection *intersection) const
{
    bool intersects = false;

    for (const Entity &entity : entities)
    {
        if (entity.intersect(ray, intersection))
            intersects = true;
    }

    return intersects;
}
