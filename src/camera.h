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
    // Screen: [0, 1]
    // Raster: [0, size in pixels]
    Transform camera_to_world;
    Transform camera_to_screen;
    Transform raster_to_screen;

    Film *film;

    Camera(const Transform &camera_to_world, float vfov, Film *film);

    float generate_ray(const CameraSample &sample, Ray *ray) const;
};
