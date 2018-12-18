#include "math.h"
#include "vector.h"
#include "point.h"
#include "ray.h"
#include "transform.h"
#include "camera.h"
#include "film.h"
#include "geometry.h"
#include "light.h"
#include "spectrum.h"

#include <fstream>
#include <sstream>
#include <memory>

static Vector3f uniform_sample_sphere(const Vector2f &u)
{
    float z = 1 - 2 * u[0];
    float r = std::sqrt(std::max((float)0, (float)1 - z * z));
    float phi = 2 * PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

static float uniform_cone_pdf(float cos_theta_max)
{
    return 1 / (2 * PI * (1 - cos_theta_max));
}

struct Entity
{
    std::shared_ptr<Geometry> geometry;
    std::shared_ptr<Light> light;

    Entity(const std::shared_ptr<Geometry> &geometry, const std::shared_ptr<Light> &light) : geometry(geometry), light(light) {}

    bool intersect(const Ray &ray, Intersection *intersection) const
    {
        if (!geometry->intersect(ray, intersection))
            return false;

        intersection->entity = this;
        return true;
    }

#if 0
    // TODO: return spectrum
    // TODO: visibility tester
    Vector3f sample(const Intersection &ref, const Vector2f &u,
                          Vector3f *wi, float *pdf) const
    {
        assert(light);

        Intersection intersection = geometry->sample(ref, u);
        *wi = normalize(intersection.p - ref.p);
        *pdf = geometry->pdf(ref, *wi);
        // TODO: visibility tester
        return light->radiance(intersection, -(*wi));
    }
#endif
};

struct Scene
{
    std::vector<Entity> entities;

    bool intersect(const Ray &ray, Intersection *intersection) const
    {
        bool intersects = false;

        for (const Entity &entity : entities)
        {
            // TODO CLEANUP: get rid of t out parameter?
            if (entity.intersect(ray, intersection))
                intersects = true;
        }

        return intersects;
    }

    void add_light(const std::shared_ptr<Geometry> geo, float emittance)
    {
        entities.push_back(Entity(geo, std::make_shared<Light>(emittance, geo)));
    }

    void add_mesh(const std::shared_ptr<Mesh> &mesh)
    {
        for (const TriangleData &t : mesh->data.tris)
            entities.push_back(Entity(std::make_shared<Triangle>(mesh, t), nullptr));
    }
};

static float uniform_float()
{
    // TODO: improve
    return (float)drand48();
}

#if 0
static float power_heuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f * f) / (f * f + g * g);
}

// TODO: spectrum
static Vector3f estimate_direct(const Light &light, const Intersection &intersection)
{
    // TODO: spectrum
    Vector3f direct(0);

    Vector2f u_light(uniform_float(), uniform_float());
    Vector3f wi;
    float light_pdf = 0;
    float scattering_pdf = 0;
    Vector3f incident = light.sample_incident_radiance(intersection, u_light, &wi, &light_pdf);
    if ((light_pdf > 0)/* && !incident.is_black()*/) // TODO: spectrum
    {
        Vector3f f(1); // TODO FIXME
        if (true/*!f.is_black()*/) // TODO FIXME
        {
            // TODO: scattering PDF
            // TODO: visibility
            float weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
            direct += Vector3f(f.x * incident.x, f.y * incident.y, f.z * incident.z) * weight / light_pdf;
//            direct += f * incident * weight / light_pdf;
        }
    }

    return direct;
}

// TODO: spectrum
static Vector3f uniform_sample_all_lights(const Scene &scene, const Intersection &intersection)
{
    Vector3f radiance(0);

    // TODO: keep separate list of lights?
    for (const Entity &entity : scene.entities)
    {
        const std::shared_ptr<Light> &light = entity.light;
        if (!light)
            continue;

        const int samples_count = 1;

        Vector3f direct(0);
        for (int i = 0; i < samples_count; ++i)
            direct += estimate_direct(*light, intersection);

        radiance += direct / (float)samples_count;
    }

    return radiance;
}
#endif

static Spectrum Li(const Scene &scene, const Ray &ray)
{
    Spectrum L(0);

    Intersection is;
    if (!scene.intersect(ray, &is))
    {
        // TODO: light emission
        return L;
    }

    // TODO FIXME

    return L;

#if 0
    const std::shared_ptr<Light> &light = intersection.entity->light;
    if (light)
        radiance += light->L(intersection, intersection.wo);

    radiance += uniform_sample_all_lights(scene, intersection);

    return radiance;
#endif
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

int main(int argc, char *argv[])
{
    Scene scene;
#if 1
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(10, 10, 4)), 0.5), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-1.25, 0, 0)), 0.1), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-3.75, 0, 0)), 0.03333), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 1.25, 0, 0)), 0.3), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 3.75, 0, 0)), 0.9), 1);

    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/floor.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate1.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate2.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate3.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate4.obj")));

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
