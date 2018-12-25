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
        entities.push_back(Entity(new Triangle(mesh, t), nullptr, bsdf));
}
