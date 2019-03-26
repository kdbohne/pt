#include "integrator.h"
#include "scene.h"
#include "intersection.h"
#include "bsdf.h"
#include "sampling.h"
#include "light.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"
#include "material.h"
#include "memory.h"

Spectrum Integrator::specular_common(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth, BxdfType type) const
{
    Point2f u(uniform_float(), uniform_float()); // TODO: better sampling?

    Vector3f wi;
    float pdf;
    BxdfType sampled_type;
    Spectrum f = its.bsdf->sample_f(its, its.wo, &wi, u, &pdf, type, &sampled_type);

    if ((pdf > 0) && !f.is_black() && (abs_dot(wi, its.n) != 0))
        return f * Li(scene, its.spawn_ray(wi), arena, depth + 1) * abs_dot(wi, its.n) / pdf;

    return Spectrum(0);
}

Spectrum Integrator::specular_reflect(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth) const
{
    BxdfType type = (BxdfType)(BSDF_SPECULAR | BSDF_REFLECTION);
    return specular_common(scene, its, ray, arena, depth, type);
}

Spectrum Integrator::specular_transmit(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth) const
{
    BxdfType type = (BxdfType)(BSDF_SPECULAR | BSDF_TRANSMISSION);
    return specular_common(scene, its, ray, arena, depth, type);
}

Spectrum WhittedIntegrator::Li(const Scene &scene, const Ray &ray, MemoryArena &arena, int depth) const
{
    Spectrum L(0);

    Intersection its;
    if (!scene.intersect(ray, &its))
    {
        // TODO: background radiance
        return Spectrum(0);
    }

    if (!its.entity->material)
        return Spectrum(0);

    its.entity->material->evaluate_surface(ray, arena, &its);

#if 1
    if (!its.bsdf)
        return Spectrum(0);
#else
    // TODO: fix this
    if (!its.entity->bsdf)
        return Li(scene, its.spawn_ray(ray.d), depth);
#endif

    // TODO: does this ever change?
    Vector3f wo = -ray.d;
    Normal3f n = its.n;

    // TODO: clean this up, add Intersection::Le() helper?
    // Hit an area light.
    if (its.entity->area_light)
        L += its.entity->area_light->L(its, wo);

    // Sum individual light contributions.
    for (const Light *light : scene.lights)
    {
        Point2f u(uniform_float(), uniform_float());

        Vector3f wi;
        float pdf;
        VisibilityTest vis;
        Spectrum Li = light->sample_Li(its, u, &wi, &pdf, &vis);

        if ((pdf == 0) || Li.is_black())
            continue;

        Spectrum f = its.bsdf->f(its, wo, wi);

        if (!f.is_black() && vis.unoccluded(scene))
            L += f * Li * abs_dot(wi, n) / pdf;
    }

    if (depth + 1 < max_depth)
    {
        L += specular_reflect(scene, its, ray, arena, depth);
        L += specular_transmit(scene, its, ray, arena, depth);
    }

    return L;
}


void WhittedIntegrator::render(const Scene &scene, const Camera &camera) const
{
    Film *film = camera.film;

    MemoryArena arena;

    for (int y = 0; y < film->resolution.y; ++y)
    {
        for (int x = 0; x < film->resolution.x; ++x)
        {
            Spectrum L(0);

            Point2i pixel(x, y);
            sampler->start_pixel(pixel);
            do
            {
                CameraSample cs = sampler->get_camera_sample(pixel);

                Ray ray;
                camera.generate_ray(cs, &ray);

                L += Li(scene, ray, arena);
            }
            while (sampler->start_next_sample());

            // TODO: proper weighting
            L /= (float)sampler->samples_per_pixel;

            film->set_pixel(x, y, L);

            arena.reset();

            float progress = (float)(y * film->resolution.x + x) / (float)(film->resolution.x * film->resolution.y);
            std::printf("\r%0.2f%%", progress * 100);
        }
    }
    std::printf("\n");

    // TODO: configurable output path/type
    camera.film->write_ppm("out.ppm");
//    camera.film->write_exr("out.exr");
}
