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
#include "parser.h"

#include <fstream>
#include <sstream>
#include <memory>

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
    Bsdf *light_bsdf = new Bsdf();
    light_bsdf->add(new LambertianReflection(Spectrum(0, 0, 0)));

    std::string path = "asset/simple.pbrt";

    Scene scene;
    Camera *camera;
    Integrator *integrator;
    if (!parse_pbrt(path, &scene, &camera, &integrator))
    {
        fatal("Failed to read input file: %s", path.c_str());
        return 1;
    }

    for (Light *light : scene.lights)
        light->preprocess(scene);

    render(scene, *camera, *integrator);
    camera->film->write_ppm("out.ppm");

    return 0;
}
