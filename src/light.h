#pragma once

#include "intersection.h"
#include "spectrum.h"
#include "vector.h"
#include "point.h"

struct Scene;
struct Transform;

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
    virtual void preprocess(const Scene &scene) {}

    virtual Spectrum Le(const Vector3f &wo) const;
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

    DirectionalLight(const Transform &light_to_world, const Spectrum &L, const Vector3f &w);

    void preprocess(const Scene &scene) override;

    Spectrum sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const override;
    Spectrum power() const override;
};
