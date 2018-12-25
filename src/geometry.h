#pragma once

#include "transform.h"

#include <vector>
#include <string>

struct Ray;

struct Geometry
{
    Transform object_to_world; // TODO: cache?

    Geometry(const Transform &object_to_world);
    virtual ~Geometry() {}

    virtual bool intersect(const Ray &ray, Intersection *intersection) const = 0;
};

struct Sphere : public Geometry
{
    float radius;

    Sphere(const Transform &object_to_world, float radius);

    bool intersect(const Ray &ray, Intersection *intersection) const override;
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

struct Mesh
{
    MeshData data;

    Mesh(const Transform &transform, const MeshData &data);
};

struct Triangle : public Geometry
{
    const Mesh *mesh;

    // TODO: reduce memory usage here
    int pi[3];
    int ni[3];
    int uvi[3];

    Triangle(const Mesh *mesh, const TriangleData &data);

    bool intersect(const Ray &ray, Intersection *intersection) const override;
};
