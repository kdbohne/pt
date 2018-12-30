#include "intersection.h"

// TODO: move this somewhere else
static constexpr float SHADOW_EPSILON = 0.0001;

Intersection::Intersection(const Point3f &p)
    : p(p), n(Normal3f()), wo(Vector3f()), time(0), entity(nullptr)
{
}

Ray Intersection::spawn_ray(const Vector3f &d) const
{
    // TODO: OffsetRayOrigin()?
    return Ray(p, d, INFINITY, time);
}

Ray Intersection::spawn_ray_to(const Intersection &p1) const
{
    // TODO: improve, two sqrt() done here
    Vector3f d = normalize(p1.p - p);
    Point3f p_biased = p + d * SHADOW_EPSILON;
    float dist = distance(p1.p, p_biased);

    return Ray(p_biased, d, dist, time);

    // NOTE: this does not work if p is not biased (OffsetRayOrigin in pbrt).
#if 0
    Vector3f d = p1.p - p;
    return Ray(p, d, 1 - SHADOW_EPSILON, time);
#endif
}
