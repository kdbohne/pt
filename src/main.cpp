#include "scene.h"
#include "camera.h"
#include "integrator.h"
#include "parser.h"
#include "light.h"

#include <cstdio>
#include <string>

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

    integrator->render(scene, *camera);

    return 0;
}
