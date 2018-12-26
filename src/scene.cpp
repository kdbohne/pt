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

void Scene::add_mesh(Mesh *mesh, Bsdf *bsdf)
{
    for (const TriangleData &t : mesh->data.tris)
    {
        Triangle *triangle = new Triangle(mesh->object_to_world, mesh->world_to_object, mesh, t);
        entities.push_back(Entity(triangle, nullptr, bsdf));
    }
}
