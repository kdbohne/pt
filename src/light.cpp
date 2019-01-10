#include "light.h"
#include "transform.h"
#include "intersection.h"
#include "math.h"
#include "scene.h"
#include "geometry.h"
#include "texture.h"
#include "sampling.h"

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

InfiniteAreaLight::InfiniteAreaLight(const Transform &light_to_world, const Spectrum &L, int samples_count, const std::string &texture_path)
    : Light(light_to_world, samples_count)
{
    Point2i resolution;
    RgbSpectrum *texels = nullptr;
    if (!texture_path.empty())
    {
        texels = read_image(texture_path, &resolution);
        for (int i = 0; i < resolution.x * resolution.y; ++i)
            texels[i] *= L;
    }
    else
    {
        resolution = Point2i(1, 1);
        texels = new RgbSpectrum[1];
        texels[0] = L;
    }

    Lmap = new Mipmap(resolution, texels);

    int width = resolution.x;
    int height = resolution.y;
    float filter = (float)1 / (float)std::max(width, height);

    float *img = new float[width * height];
    for (int v = 0; v < height; ++v)
    {
        float vp = (float)v / (float)height;
        float sin_theta = std::sin(PI * (v + 0.5) / (float)height);
        for (int u = 0; u < width; ++u)
        {
            float up = (float)u / (float)width;
            img[v * width + u] = Lmap->lookup(Point2f(up, vp), filter).y() * sin_theta;
        }
    }

    distribution = new Distribution2d(img, width, height);
}

void InfiniteAreaLight::preprocess(const Scene &scene)
{
    // TODO CLEANUP: this is duplicated from DirectionalLight::preprocess()
    Bounds3f bounds;
    for (const Entity &e : scene.entities)
        bounds = Union(e.geometry->world_bounds(), bounds);

    bounds.bounding_sphere(&world_center, &world_radius);
}

Spectrum InfiniteAreaLight::Le(const Ray &ray) const
{
    Transform world_to_light = inverse(light_to_world); // TODO: store in Light

    Vector3f w = normalize(world_to_light * ray.d);
    Point2f st(spherical_phi(w) * INV_2_PI,
               spherical_theta(w) * INV_2_PI);

    return Spectrum(Lmap->lookup(st), SpectrumType::ILLUMINANT);
}

Spectrum InfiniteAreaLight::sample_Li(const Intersection &ref, const Point2f &u, Vector3f *wi, float *pdf, VisibilityTest *vis) const
{
    float map_pdf;
    Point2f uv = distribution->sample_continuous(u, &map_pdf);
    if (map_pdf == 0)
        return Spectrum(0);

    float theta = uv[1] * PI;
    float phi = uv[0] * 2 * PI;

    float cos_theta = std::cos(theta);
    float sin_theta = std::sin(theta);
    float cos_phi = std::cos(phi);
    float sin_phi = std::sin(phi);

    *wi = light_to_world * Vector3f(sin_theta * cos_phi,
                                    sin_theta * sin_phi,
                                    cos_theta);

    *pdf = map_pdf / (2 * PI * PI * sin_theta);
    if (sin_theta == 0)
        *pdf = 0;

    // TODO: set p1's time to ref.time
    *vis = VisibilityTest(ref, Intersection(ref.p + *wi * (2 * world_radius)));

    return Spectrum(Lmap->lookup(uv), SpectrumType::ILLUMINANT);
}

Spectrum InfiniteAreaLight::power() const
{
    Spectrum L = Spectrum(Lmap->lookup(Point2f(0.5, 0.5), 0.5), SpectrumType::ILLUMINANT);
    return PI * world_radius * world_radius * L;
}
