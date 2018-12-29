#include "scene.h"
#include "camera.h"
#include "integrator.h"
#include "film.h"
#include "ray.h"
#include "spectrum.h"
#include "parser.h"
#include "light.h"

static void render(const Scene &scene, const Camera &camera, const Integrator &integrator)
{
    Film *film = camera.film;

    for (int y = 0; y < film->resolution.y; ++y)
    {
        for (int x = 0; x < film->resolution.x; ++x)
        {
            CameraSample cs;
            cs.film_pos.x = (float)x / (float)(film->resolution.x - 1);
            cs.film_pos.y = 1 - (float)y / (float)(film->resolution.y - 1);
            cs.time = 0; // TODO

            Ray ray;
            camera.generate_ray(cs, &ray);

            Spectrum L = integrator.Li(scene, ray);
            film->set_pixel(x, y, L);
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::printf("Usage: ./pt <scene>\n");
        return 0;
    }

    std::string path = argv[1];

    Scene scene;
    Camera *camera;
    Integrator *integrator;
    if (!parse_pbrt(path, &scene, &camera, &integrator))
        return 1;

    for (Light *light : scene.lights)
        light->preprocess(scene);

    render(scene, *camera, *integrator);
    camera->film->write_ppm("out.ppm");

    return 0;
}
