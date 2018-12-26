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
#include "integrator.h"

#include <fstream>
#include <sstream>
#include <memory>

static void render(const Scene &scene, const Camera &camera, const Integrator &integrator, Film &film)
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

            Spectrum L = integrator.Li(scene, ray);
            film.set_pixel(x, y, L);
        }
    }
}

static Bsdf *make_matte_material(const Spectrum &Kd)
{
    Bsdf *bsdf = new Bsdf();

    if (!Kd.is_black())
        bsdf->add(new LambertianReflection(Kd));

    return bsdf;
}

static Bsdf *make_plastic_material(const Spectrum &Kd, const Spectrum &Ks, float roughness)
{
    Bsdf *bsdf = new Bsdf();

    if (!Kd.is_black())
        bsdf->add(new LambertianReflection(Kd));

    if (!Ks.is_black())
    {
        MicrofacetDistribution *distribution = new TrowbridgeReitzDistribution(roughness, roughness);
        Fresnel *fresnel = new FresnelDielectric(1, 1.5);
        MicrofacetReflection *specular = new MicrofacetReflection(Ks, distribution, fresnel);
        bsdf->add(specular);
    }

    return bsdf;
}

int main(int argc, char *argv[])
{
    Bsdf *light_bsdf = new Bsdf();
    light_bsdf->add(new LambertianReflection(Spectrum(0, 0, 0)));

    Scene scene;
    scene.lights.push_back(new DirectionalLight(Transform(), Spectrum(10, 10, 10), normalize(Point3f(0, 10, 10) - Point3f(0, 0, -10))));
#if 1
    scene.lights.push_back(new PointLight(translate(Vector3f(   10, 10, 4)), Spectrum(    800,     800,     800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(-1.25,  0, 0)), Spectrum(    100,     100,     100)));
    scene.lights.push_back(new PointLight(translate(Vector3f(-3.75,  0, 0)), Spectrum(901.803, 901.803, 901.803)));
    scene.lights.push_back(new PointLight(translate(Vector3f( 1.25,  0, 0)), Spectrum(11.1111, 11.1111, 11.1111)));
    scene.lights.push_back(new PointLight(translate(Vector3f( 3.75,  0, 0)), Spectrum(1.23457, 1.23457, 1.23457)));
#endif

    scene.add_mesh(new Mesh(Transform(), Transform(), load_obj("asset/veach_mi/floor.obj")), make_matte_material(Spectrum(0.4, 0.4, 0.4)));
    scene.add_mesh(new Mesh(Transform(), Transform(), load_obj("asset/veach_mi/plate1.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.005));
    scene.add_mesh(new Mesh(Transform(), Transform(), load_obj("asset/veach_mi/plate2.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.02));
    scene.add_mesh(new Mesh(Transform(), Transform(), load_obj("asset/veach_mi/plate3.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.05));
    scene.add_mesh(new Mesh(Transform(), Transform(), load_obj("asset/veach_mi/plate4.obj")), make_plastic_material(Spectrum(0.07, 0.09, 0.13), Spectrum(1, 1, 1), 0.1));

    Film film(Point2i(1152, 864) / 2);
    Camera camera(look_at(Vector3f(0, 2, 15), Vector3f(0, 1.69521, 14.0476), Vector3f(0, 0.952421, -0.304787)), radians(28), &film);
    WhittedIntegrator integrator(5);

#if 0
    Film film(Point2i(640, 360));
    Camera camera(look_at(Vector3f(8, 2, 3), Vector3f(0, 0, 0), Vector3f(0, 1, 0)), radians(40), &film);
    WhittedIntegrator integrator(5);

    scene.lights.push_back(new PointLight(translate(Vector3f(0, 10, 0)), Spectrum(800, 800, 800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(-5, 8, 3)), Spectrum(800, 800, 800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(5, 14, 0)), Spectrum(800, 800, 800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(2, 18, -3)), Spectrum(800, 800, 800)));
    scene.lights.push_back(new PointLight(translate(Vector3f(0, 100, 0)), Spectrum(80000, 80000, 80000)));

    scene.entities.push_back(Entity(new Sphere(translate(Vector3f(0, 0, 0)), 1000), nullptr, make_matte_material(Spectrum(0.4, 0.4, 0.4))));
    {
        scene.entities.push_back(Entity(new Sphere(translate(Vector3f(0, -1000, 0)), 1000), nullptr, make_matte_material(Spectrum(0.4, 0.4, 0.4))));

        for (int a = -11; a < 11; ++a)
        {
            for (int b = -11; b < 11; ++b)
            {
                float choice = drand48();
                Vector3f center(a + 0.9 * drand48(), 0.2, b + 0.9 * drand48());
                if ((center - Vector3f(4.0, 0.2, 0.0)).length() > 0.9)
                {
                    if (choice < 0.8)
                        scene.entities.push_back(Entity(new Sphere(translate(center), 0.2), nullptr, make_matte_material(Spectrum(drand48() * drand48(), drand48() * drand48(), drand48() * drand48()))));
                    else if (choice < 0.95)
                        scene.entities.push_back(Entity(new Sphere(translate(center), 0.2), nullptr, make_matte_material(Spectrum(0.5 * (1.0 + drand48()), 0.5 * (1.0 + drand48()), 0.5 * (1.0 + drand48())))));
                    else
                        scene.entities.push_back(Entity(new Sphere(translate(center), 0.2), nullptr, make_plastic_material(Spectrum(0.5, 0.5, 0.5), Spectrum(1, 1, 1), 0.1)));
                }
            }
        }

        scene.entities.push_back(Entity(new Sphere(translate(Vector3f( 0, 1, 0)), 1), nullptr, make_plastic_material(Spectrum(0.2, 0.2, 0.3), Spectrum(1, 1, 1), 0.05)));
        scene.entities.push_back(Entity(new Sphere(translate(Vector3f(-4, 1, 0)), 1), nullptr, make_matte_material(Spectrum(0.4, 0.2, 0.1))));
        scene.entities.push_back(Entity(new Sphere(translate(Vector3f( 4, 1, 0)), 1), nullptr, make_matte_material(Spectrum(0.7, 0.6, 0.5))));
    }
#endif

    for (Light *light : scene.lights)
        light->preprocess(scene);

    render(scene, camera, integrator, film);
    film.write_ppm("out.ppm");

    return 0;
}
