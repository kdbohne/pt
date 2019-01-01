#include "integrator.h"
#include "scene.h"
#include "intersection.h"
#include "bsdf.h"
#include "sampling.h"
#include "light.h"
#include "camera.h"
#include "film.h"
#include "sampler.h"

Spectrum Integrator::specular_common(const Scene &scene, const Intersection &its, const Ray &ray, int depth, BxdfType type) const
{
    Point2f u(uniform_float(), uniform_float()); // TODO: better sampling?

    Vector3f wi;
    float pdf;
    BxdfType sampled_type;
    Spectrum f = its.entity->bsdf->sample_f(its, its.wo, &wi, u, &pdf, type, &sampled_type);

    if ((pdf > 0) && !f.is_black() && (abs_dot(wi, its.n) != 0))
        return f * Li(scene, its.spawn_ray(wi), depth + 1) * abs_dot(wi, its.n) / pdf;

    return Spectrum(0);
}

Spectrum Integrator::specular_reflect(const Scene &scene, const Intersection &its, const Ray &ray, int depth) const
{
    BxdfType type = (BxdfType)(BSDF_SPECULAR | BSDF_REFLECTION);
    return specular_common(scene, its, ray, depth, type);
}

Spectrum Integrator::specular_transmit(const Scene &scene, const Intersection &its, const Ray &ray, int depth) const
{
    BxdfType type = (BxdfType)(BSDF_SPECULAR | BSDF_TRANSMISSION);
    return specular_common(scene, its, ray, depth, type);
}

Spectrum WhittedIntegrator::Li(const Scene &scene, const Ray &ray, int depth) const
{
    Spectrum L(0);

    Intersection its;
    if (!scene.intersect(ray, &its))
    {
        // TODO: background radiance
        return Spectrum(0);
    }

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

        Spectrum f = its.entity->bsdf->f(its, wo, wi);

        if (!f.is_black() && vis.unoccluded(scene))
            L += f * Li * abs_dot(wi, n) / pdf;
    }

    if (depth + 1 < max_depth)
    {
        L += specular_reflect(scene, its, ray, depth);
        L += specular_transmit(scene, its, ray, depth);
    }

    return L;
}


void WhittedIntegrator::render(const Scene &scene, const Camera &camera) const
{
    Film *film = camera.film;

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

                L += Li(scene, ray);
            }
            while (sampler->start_next_sample());

            // TODO: proper weighting
            L /= (float)sampler->samples_per_pixel;

            film->set_pixel(x, y, L);
        }
    }

    // TODO: configurable output path/type
    camera.film->write_ppm("out.ppm");
}

#if 0
static float power_heuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f * f) / (f * f + g * g);
}

static Spectrum estimate_direct(const Scene &scene, const Intersection &its, const Light &light, const Point2f &u_scattering, bool specular = false)
{
    const std::shared_ptr<Bsdf> &bsdf = its.entity->bsdf;
    if (!bsdf)
    {
        // TODO: default material
        return Spectrum(0);
    }

    BxdfType bsdf_flags = specular ? BSDF_ALL : (BxdfType)(BSDF_ALL & ~BSDF_SPECULAR);

    Spectrum Ld(0);

    Point2f u_light(uniform_float(), uniform_float());
    Vector3f wi;
    float light_pdf = 0, scattering_pdf = 0;
    Spectrum Li = light.sample_Li(its, u_light, &wi, &light_pdf);

    if ((light_pdf > 0) && !Li.is_black())
    {
        // TODO: assuming always surface intersection; handle participating media
        Spectrum f = bsdf->f(its, its.wo, wi, bsdf_flags) * abs_dot(wi, its.n); // TODO: normal map (its.shading.n)
        scattering_pdf = bsdf->pdf(its, its.wo, wi, bsdf_flags);

        if (!f.is_black())
        {
            // TODO: visibility
            /*
            if (!visible)
                Li = Spectrum(0);
            */

            if (!Li.is_black())
            {
                // TODO: delta light
                float weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
                Ld += f * Li * weight / light_pdf;
            }
        }
    }

#if 0
//    if (!is_delta_light(light.flags)) // TODO FIXME
    {
        // TODO: assuming always surface intersection; handle participating media
        BxdfType sampled_type;
        Spectrum f = bsdf->sample_f(its, its.wo, &wi, u_scattering, &scattering_pdf, bsdf_flags, &sampled_type) * abs_dot(wi, its.n); // TODO: normal map (its.shading.n)

        bool sampled_specular = sampled_type & BSDF_SPECULAR;

        if (!f.is_black() && (scattering_pdf > 0))
        {
            float weight = 1;
            if (!sampled_specular)
            {
                light_pdf = light.pdf_Li(its, wi);
                if (light_pdf == 0)
                    return Ld;

                weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
            }

            Ray ray = its.spawn_ray(wi);

            Spectrum Tr(1); // NOTE: this would be set by participating media intersect_Tr()
            Intersection light_is;
            if (scene.intersect(ray, &light_is))
            {
                Li = Spectrum(0);
                if (light_is.entity->light) // TODO: check if area light
                    Li = light_is.Le(-wi);
            }
            else
            {
                Li = light.Le(ray);
            }

            if (!Li.is_black())
                Ld += f * Li * Tr * weight / scattering_pdf;
        }
    }
#endif

    return Ld;
}

// TODO: spectrum
static Spectrum uniform_sample_all_lights(const Scene &scene, const Intersection &intersection)
{
    Spectrum L(0);

    // TODO: keep separate list of lights?
    for (const Entity &entity : scene.entities)
    {
        const std::shared_ptr<Light> &light = entity.light;
        if (!light)
            continue;

        Point3f p_light = entity.geometry->object_to_world * Point3f(0, 0, 0);
        float d = distance(intersection.p, p_light);

        const int samples_count = 1; // TODO: ?

        Spectrum Ld(0);
        for (int i = 0; i < samples_count; ++i)
        {
            Point2f u_scattering(uniform_float(), uniform_float()); // TODO: better sampling?
            Ld += estimate_direct(scene, intersection, *light, u_scattering);
        }

        L += Ld / (float)samples_count;
    }

    return L;
}
#endif

#if 0
static Spectrum Li(const Scene &scene, const Ray &ray)
{
    Spectrum L(0);

    Intersection its;
    if (!scene.intersect(ray, &its))
    {
        // TODO: light emission
        return L;
    }

#if 0
    if (!its.bsdf)
        return Li(scene, its.spawn_ray(ray.d));
#endif

    const std::shared_ptr<Light> &light = its.entity->light;
    if (light)
        L += light->L(its, its.wo);

    L += uniform_sample_all_lights(scene, its);

    return L;
}
#endif
