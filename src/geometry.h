#pragma once

#include "transform.h"
#include "bounds.h"

#include <vector>
#include <string>

struct Ray;

struct Geometry
{
    Transform object_to_world, world_to_object; // TODO: cache?

    Geometry(const Transform &object_to_world, const Transform &world_to_object);
    virtual ~Geometry() {}

    virtual Bounds3f object_bounds() const = 0;
    virtual Bounds3f world_bounds() const;

    virtual float area() const = 0;

    virtual bool intersect(const Ray &ray, Intersection *intersection) const = 0;

    virtual Intersection sample(const Point2f &u, float *pdf) const = 0;
    virtual Intersection sample(const Intersection &ref, const Point2f &u, float *pdf) const;
};

struct Sphere : public Geometry
{
    float radius;

    Sphere(const Transform &object_to_world, const Transform &world_to_object, float radius);

    Bounds3f object_bounds() const override;

    float area() const override;

    bool intersect(const Ray &ray, Intersection *intersection) const override;

    Intersection sample(const Point2f &u, float *pdf) const override;
    Intersection sample(const Intersection &ref, const Point2f &u, float *pdf) const override;
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
    std::vector<Normal3f> n;
    std::vector<Point2f> uv;
};

MeshData load_obj(const std::string &path);
MeshData load_ply(const std::string &path);

struct Mesh
{
    // TODO CLEANUP: this is sort of duplicated from Geometry... not sure how
    // best to organize mesh/triangle types
    Transform object_to_world, world_to_object; // TODO: cache?

    MeshData data;

    Mesh(const Transform &object_to_world, const Transform &world_to_object, const MeshData &data);
};

struct Triangle : public Geometry
{
    Mesh *mesh;

    // TODO: reduce memory usage here
    int pi[3];
    int ni[3];
    int uvi[3];

    Triangle(const Transform &object_to_world, const Transform &world_to_object, Mesh *mesh, const TriangleData &data);

    Bounds3f object_bounds() const override;
    Bounds3f world_bounds() const override;

    float area() const override;

    bool intersect(const Ray &ray, Intersection *intersection) const override;

    Intersection sample(const Point2f &u, float *pdf) const override;
};
