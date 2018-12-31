#pragma once

#include "vector.h"
#include "point.h"
#include "transform.h"
#include "film.h"

struct CameraSample
{
    Point2f p_film;
    float time;
};

struct Camera
{
    Transform camera_to_world;
    Transform camera_to_screen;
    Transform raster_to_screen;

    Film *film;

    Camera(const Transform &camera_to_world, float vfov, Film *film)
        : camera_to_world(camera_to_world),
          camera_to_screen(perspective(vfov, (float)film->resolution.x / (float)film->resolution.y, 1e-2f, 1000.0f)),
          raster_to_screen(translate(Vector3f(-1, 1, 0)) * scale((float)2 / (float)film->resolution.x, -(float)2 / (float)film->resolution.y, 1)),
          film(film)
    {
    }

    float generate_ray(const CameraSample &sample, Ray *ray) const;
};
