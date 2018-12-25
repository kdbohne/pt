#include "entity.h"
#include "intersection.h"
#include "geometry.h"

bool Entity::intersect(const Ray &ray, Intersection *intersection) const
{
    if (!geometry->intersect(ray, intersection))
        return false;

    intersection->entity = this;
    return true;
}
