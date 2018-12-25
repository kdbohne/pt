#include "math.h"
#include "vector.h"
#include "point.h"
#include "ray.h"
#include "entity.h"
#include "transform.h"
#include "camera.h"
#include "film.h"
#include "geometry.h"
#include "light.h"
#include "spectrum.h"
#include "bsdf.h"
#include "scene.h"
#include "sampling.h"

#include <fstream>
#include <sstream>
#include <memory>

#if 0
static float power_heuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f * f) / (f * f + g * g);
}

static Spectrum estimate_direct(const Scene &scene, const Intersection &is, const Light &light, const Point2f &u_scattering, bool specular = false)
{
    const std::shared_ptr<Bsdf> &bsdf = is.entity->bsdf;
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
    Spectrum Li = light.sample_Li(is, u_light, &wi, &light_pdf);

    if ((light_pdf > 0) && !Li.is_black())
    {
        // TODO: assuming always surface intersection; handle participating media
        Spectrum f = bsdf->f(is, is.wo, wi, bsdf_flags) * abs_dot(wi, is.n); // TODO: normal map (is.shading.n)
        scattering_pdf = bsdf->pdf(is, is.wo, wi, bsdf_flags);

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
        Spectrum f = bsdf->sample_f(is, is.wo, &wi, u_scattering, &scattering_pdf, bsdf_flags, &sampled_type) * abs_dot(wi, is.n); // TODO: normal map (is.shading.n)

        bool sampled_specular = sampled_type & BSDF_SPECULAR;

        if (!f.is_black() && (scattering_pdf > 0))
        {
            float weight = 1;
            if (!sampled_specular)
            {
                light_pdf = light.pdf_Li(is, wi);
                if (light_pdf == 0)
                    return Ld;

                weight = power_heuristic(1, scattering_pdf, 1, light_pdf);
            }

            Ray ray = is.spawn_ray(wi);

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

    Intersection is;
    if (!scene.intersect(ray, &is))
    {
        // TODO: light emission
        return L;
    }

#if 0
    if (!is.bsdf)
        return Li(scene, is.spawn_ray(ray.d));
#endif

    const std::shared_ptr<Light> &light = is.entity->light;
    if (light)
        L += light->L(is, is.wo);

    L += uniform_sample_all_lights(scene, is);

    return L;
}
#endif

#if 0
static Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0);

static Spectrum specular_reflect(const Scene &scene, const Intersection &is, const Ray &ray, int depth)
{
    BxdfType type = (BxdfType)(BSDF_REFLECTION | BSDF_SPECULAR);
    Point2f u(uniform_float(), uniform_float()); // TODO: better sampling?

    Vector3f wi;
    float pdf;
    Spectrum f = is.entity->bsdf->sample_f(is, is.wo, &wi, u, &pdf, type, nullptr);

    if ((pdf > 0) && !f.is_black() && (abs_dot(wi, is.n) != 0))
        return f * Li(scene, ray, depth + 1) * abs_dot(wi, is.n) / pdf;

    return Spectrum(0);
}

static Spectrum specular_transmit(const Scene &scene, const Intersection &is, const Ray &ray, int depth)
{
    // TODO FIXME
    return Spectrum(0);
}

static Spectrum Li(const Scene &scene, const Ray &ray, int depth)
{
    Spectrum L(0);

    Intersection is;
    if (!scene.intersect(ray, &is))
    {
        // TODO: separate lights list?
        for (const Entity &entity : scene.entities)
        {
            const std::shared_ptr<Light> &light = entity.light;
            if (!light)
                continue;

            L += light->Le(ray);
        }

        return L;
    }

    if (!is.entity->bsdf)
        return Li(scene, is.spawn_ray(ray.d), depth);

    L += is.Le(is.wo);

    for (const Entity &entity : scene.entities)
    {
        const std::shared_ptr<Light> &light = entity.light;
        if (!light)
            continue;

        Point2f u(uniform_float(), uniform_float()); // TODO: better sampling?
        Vector3f wi;
        float pdf;
        Spectrum Li = light->sample_Li(is, u, &wi, &pdf);
        if (Li.is_black() || (pdf == 0))
            continue;

        Spectrum f = is.entity->bsdf->f(is, is.wo, wi);
        if (!f.is_black() /*&& !occluded*/)
            L += f * Li * abs_dot(wi, is.n) / pdf;
    }

    static const int MAX_DEPTH = 5;
    if (depth + 1 < MAX_DEPTH)
    {
        L += specular_reflect(scene, is, ray, depth);
        L += specular_transmit(scene, is, ray, depth);
    }

    return L;
}
#endif

static Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0)
{
    Spectrum L(0);

    Intersection is;
    if (!scene.intersect(ray, &is))
    {
        // TODO: background radiance
        return Spectrum(0);
    }

    // TODO: does this ever change?
    Vector3f wo = -ray.d;
    Normal3f n = is.n;

    // Hit an area light.
    if (is.entity->light)
        L += is.entity->light->Le(wo);

    // Sum individual light contributions.
    for (const Light *light : scene.lights)
    {
        Point2f u(uniform_float(), uniform_float());

        Vector3f wi;
        float pdf;
        VisibilityTest vis;
        Spectrum Li = light->sample_Li(is, u, &wi, &pdf, &vis);

        if ((pdf == 0) || Li.is_black())
            continue;

        Spectrum f = is.entity->bsdf->f();

        if (!f.is_black() && vis.unoccluded(scene))
            L += f * Li * abs_dot(wi, n) / pdf;
    }

    return L;
}

static void render(const Scene &scene, const Camera &camera, Film &film)
{
    for (int y = 0; y < film.resolution.y; ++y)
    {
        for (int x = 0; x < film.resolution.x; ++x)
        {
            CameraSample cs;
            cs.film_pos.x = (float)x / (float)(film.resolution.x - 1);
            cs.film_pos.y = 1 - (float)y / (float)(film.resolution.y - 1);
            cs.time = 0; // TODO

            Ray ray;
            camera.generate_ray(cs, &ray);

            Spectrum L = Li(scene, ray);
            film.set_pixel(x, y, L);
        }
    }
}

static Bsdf *make_matte_material(const Spectrum &Kd)
{
    Bsdf *bsdf = new Bsdf;

    if (!Kd.is_black())
        bsdf->add(new LambertianReflection(Kd));

    return bsdf;
}

#if 0
static std::shared_ptr<Bsdf> make_plastic_material(const Spectrum &Kd, const Spectrum &Ks, float roughness)
{
    std::shared_ptr<Bsdf> bsdf = std::make_shared<Bsdf>();

    if (!Kd.is_black())
        bsdf->add(std::make_shared<LambertianReflection>(Kd));

    // TODO FIXME
#if 0
    if (!Ks.is_black())
    {
        std::shared_ptr<MicrofacetDistribution> distribution = std::make_shared<TrowbridgeReitzDistribution>(roughness, roughness);
        std::shared_ptr<Fresnel> fresnel = std::make_shared<FresnelDielectric>(1, 1.5);
        std::shared_ptr<MicrofacetReflection> specular = std::make_shared<MicrofacetReflection>(Ks, distribution, fresnel);
        bsdf->add(specular);
    }
#endif

    return bsdf;
}
#endif

int main(int argc, char *argv[])
{
    Bsdf *light_bsdf = new Bsdf();
    light_bsdf->add(new LambertianReflection(Spectrum(0, 0, 0)));

    Scene scene;
    scene.lights.push_back(new PointLight(translate(Vector3f(   10, 10, 4)), Spectrum(    800,     800,     800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(-1.25,  0, 0)), Spectrum(    100,     100,     100)));
    scene.lights.push_back(new PointLight(translate(Vector3f(-3.75,  0, 0)), Spectrum(901.803, 901.803, 901.803)));
    scene.lights.push_back(new PointLight(translate(Vector3f( 1.25,  0, 0)), Spectrum(11.1111, 11.1111, 11.1111)));
    scene.lights.push_back(new PointLight(translate(Vector3f( 3.75,  0, 0)), Spectrum(1.23457, 1.23457, 1.23457)));

    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/floor.obj")), make_matte_material(Spectrum(0.4, 0.4, 0.4)));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate1.obj")), make_matte_material(Spectrum(0.07, 0.09, 0.13)));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate2.obj")), make_matte_material(Spectrum(0.07, 0.09, 0.13)));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate3.obj")), make_matte_material(Spectrum(0.07, 0.09, 0.13)));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate4.obj")), make_matte_material(Spectrum(0.07, 0.09, 0.13)));
#if 0
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate1.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.005));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate2.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.02));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate3.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.05));
    scene.add_mesh(new Mesh(Transform(), load_obj("asset/veach_mi/plate4.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.1));
#endif

    Film film(Point2i(1152, 864) / 2);
    Camera camera(look_at(Vector3f(0, 2, 15), Vector3f(0, 1.69521, 14.0476), Vector3f(0, 0.952421, -0.304787)), radians(28), &film);

#if 0
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(10, 10, 4)), 0.5), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-1.25, 0, 0)), 0.1), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-3.75, 0, 0)), 0.03333), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 1.25, 0, 0)), 0.3), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 3.75, 0, 0)), 0.9), 1);

    std::shared_ptr<Bsdf> bsdf = std::make_shared<Bsdf>();
    bsdf->add(std::make_shared<LambertianReflection>(Spectrum(0.5)));

    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/floor.obj")), bsdf);
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate1.obj")), bsdf);
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate2.obj")), bsdf);
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate3.obj")), bsdf);
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate4.obj")), bsdf);

    Film film(Point2i(768, 512) / 2);
    Camera camera(look_at(Vector3f(0, 2, 15), Vector3f(0, -2, 2.5), Vector3f(0, 1, 0)), radians(28), &film);
#endif

#if 0
    Film film(Point2i(640, 360));
    Camera camera(look_at(Vector3f(8, 2, 3), Vector3f(0, 0, 0), Vector3f(0, 1, 0)), radians(40), &film);

    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(0, 0, 0)), 1), nullptr));
    {
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(0, -1000, 0)), 1000), nullptr));

        for (int a = -11; a < 11; ++a)
        {
            for (int b = -11; b < 11; ++b)
            {
                float choice = drand48();
                Vector3f center(a + 0.9 * drand48(), 0.2, b + 0.9 * drand48());
                if ((center - Vector3f(4.0, 0.2, 0.0)).length() > 0.9)
                {
                    if (choice < 0.8)
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                    else if (choice < 0.95)
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                    else
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                }
            }
        }

        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 0, 1, 0)), 1), nullptr));
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(-4, 1, 0)), 1), nullptr));
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 4, 1, 0)), 1), nullptr));
    }
#endif

    render(scene, camera, film);
    film.write_ppm("out.ppm");

    return 0;
}
