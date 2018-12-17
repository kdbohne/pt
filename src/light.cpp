#include "light.h"
#include "geometry.h"

Vector3f Light::L(const Intersection &is, const Vector3f &w) const
{
    return (dot(is.n, w) > 0) ? Lemit : Vector3f(0);
}

#if 0
Vector3f Light::sample_Li(const Intersection &ref, const Vector2f &u, Vector3f *wi, float *pdf) const
{
    Intersection geo_is = geometry->sample(ref, u, pdf);
    *wi = normalize(geo_is.p - ref.p);
    // TODO: visibility tester
    return L(geo_inter, -(*wi));
}
#endif
