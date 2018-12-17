#pragma once

#include "transform.h"

#include <memory>
#include <vector>

struct Mesh;

struct Geometry
{
    Transform object_to_world; // TODO: cache?

    Geometry(const Transform &object_to_world) : object_to_world(object_to_world) {}
    virtual ~Geometry() {}

    virtual bool intersect(const Ray &ray, Intersection *intersection) const = 0;

#if 0
    virtual Intersection sample(const Vector2f &u, float *pdf) const = 0;
    virtual Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const = 0;

    virtual float area() const = 0;
    virtual float pdf(const Intersection &ref, const Vector3f &wi) const;
#endif
};

struct Sphere : public Geometry
{
    float radius;

    Sphere();
    Sphere(const Transform &object_to_world, float radius);

    bool intersect(const Ray &ray, Intersection *intersection) const override;

#if 0
    Intersection sample(const Vector2f &u, float *pdf) const override;
    Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const override;

    float area() const override;
    float pdf(const Intersection &ref, const Vector3f &wi) const override;
#endif
};

struct TriangleData
{
    int pi[3];
    int ni[3];
    int uvi[3];
};

struct MeshData
{
    std::vector<TriangleData> tris;
    std::vector<Point3f> p;
    std::vector<Vector3f> n;
};

MeshData load_obj(const std::string &path);

struct Mesh
{
    MeshData data;

    Mesh(const Transform &transform, const MeshData &data);
};

struct Triangle : public Geometry
{
    std::shared_ptr<Mesh> mesh;
    // TODO: reduce memory usage here
    int pi[3];
    int ni[3];

    Triangle(const std::shared_ptr<Mesh> &mesh, const TriangleData &data);

    bool intersect(const Ray &ray, Intersection *intersection) const override;

#if 0
    Intersection sample(const Vector2f &u, float *pdf) const override;
    Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const override;

    float area() const override;
    float pdf(const Intersection &ref, const Vector3f &wi) const override;
#endif
};
