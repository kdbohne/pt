#pragma once

#include "intersection.h"
#include "spectrum.h"
#include "vector.h"
#include "point.h"
#include "transform.h"
#include "geometry.h"

struct Scene;
template<typename T> struct Texture;
//template<typename T> struct Mipmap;
struct Distribution2d;

class VisibilityTest
{
public:
    VisibilityTest() {}
    VisibilityTest(const Intersection &p0, const Intersection &p1);

    bool unoccluded(const Scene &scene);

private:
    Intersection p0, p1;
};

struct Light
{
    Transform light_to_world;
    int samples_count;

    Light(Transform light_to_world, int samples_count = 1)
        : light_to_world(light_to_world), samples_count(samples_count) {}

    virtual void preprocess(const Scene &scene) {}

    virtual Spectrum Le(const Ray &ray) const;
    virtual Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const = 0;
    virtual Spectrum power() const = 0;
};

struct PointLight : public Light
{
    Point3f p;
    Spectrum I;

    PointLight(const Transform &light_to_world, const Spectrum &I);

    Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const override;
    Spectrum power() const override;
};

struct DirectionalLight : public Light
{
    Spectrum L;
    Vector3f w_light;

    Point3f world_center;
    float world_radius;

    DirectionalLight(const Transform &light_to_world, const Spectrum &L, const Vector3f &w)
        : Light(light_to_world), L(L), w_light(normalize(light_to_world * w)) {}

    void preprocess(const Scene &scene) override;

    Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const override;
    Spectrum power() const override;
};

struct AreaLight : public Light
{
    AreaLight(const Transform &light_to_world, int samples_count) : Light(light_to_world, samples_count) {}

    virtual Spectrum L(const Intersection &its, const Vector3f &w) const = 0;
};

struct DiffuseAreaLight : public AreaLight
{
    Spectrum Lemit;
    Geometry *geometry; // TODO: move to AreaLight?
    float area;

    DiffuseAreaLight(const Transform &light_to_world, int samples_count, const Spectrum &Lemit, Geometry *geometry)
        : AreaLight(light_to_world, samples_count), Lemit(Lemit), geometry(geometry), area(geometry->area()) {}

    Spectrum L(const Intersection &its, const Vector3f &w) const override;

    Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const override;
    Spectrum power() const override;
};

// NOTE: although InfiniteAreaLight is an area light, it only inherits from
// Light rather than Light's subclass AreaLight.
struct InfiniteAreaLight : public Light
{
    Texture<RgbSpectrum> *Lmap;
//    Mipmap<RgbSpectrum> *Lmap;

    Distribution2d *distribution;

    Point3f world_center;
    float world_radius;

    InfiniteAreaLight(const Transform &light_to_world, const Spectrum &L, int samples_count, const std::string &texture_path);

    void preprocess(const Scene &scene) override;

    Spectrum Le(const Ray &ray) const override;
    Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const override;
    Spectrum power() const override;
};
