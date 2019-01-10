#include "light.h"
#include "transform.h"
#include "intersection.h"
#include "math.h"
#include "scene.h"
#include "geometry.h"

VisibilityTest::VisibilityTest(const Intersection &p0, const Intersection &p1)
    : p0(p0), p1(p1)
{
}

bool VisibilityTest::unoccluded(const Scene &scene)
{
    // TODO: intersect_p()
    Intersection is;
    return !scene.intersect(p0.spawn_ray_to(p1), &is);
}

Spectrum Light::Le(const Ray &ray) const
{
    return Spectrum(0);
}

PointLight::PointLight(const Transform &light_to_world, const Spectrum &I)
    : Light(light_to_world), p(light_to_world * Point3f(0, 0, 0)), I(I)
{
}

Spectrum PointLight::sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const
{
    *wi = normalize(p - ref.p);
    *pdf = 1;

    *vis = VisibilityTest(ref, Intersection(p));

    return I / distance_squared(p, ref.p);
}

Spectrum PointLight::power() const
{
    return 4 * PI * I;
}

void DirectionalLight::preprocess(const Scene &scene)
{
    Bounds3f bounds;
    for (const Entity &e : scene.entities)
        bounds = Union(e.geometry->world_bounds(), bounds);

    bounds.bounding_sphere(&world_center, &world_radius);
}

Spectrum DirectionalLight::sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const
{
    *wi = w_light;
    *pdf = 1;

    Point3f p_outside = ref.p + w_light * (2 * world_radius);

    // TODO: set p1's time to ref.time
    *vis = VisibilityTest(ref, Intersection(p_outside));

    return L;
}

Spectrum DirectionalLight::power() const
{
    return L * PI * world_radius * world_radius;
}

Spectrum DiffuseAreaLight::L(const Intersection &its, const Vector3f &w) const
{
    // TODO: two-sided option
    return (dot(its.n, w) > 0) ? Lemit : Spectrum(0);
}

Spectrum DiffuseAreaLight::sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const
{
    Intersection p_shape = geometry->sample(ref, u, pdf);
    if ((*pdf == 0) || ((p_shape.p - ref.p).length_squared() == 0))
    {
        *pdf = 0;
        return Spectrum(0);
    }

    *wi = normalize(p_shape.p - ref.p);
    *vis = VisibilityTest(ref, p_shape);

    return L(p_shape, -(*wi));
}

Spectrum DiffuseAreaLight::power() const
{
    // TODO: two-sided option
    return Lemit * area * PI;
}
