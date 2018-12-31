#include "camera.h"

#include <stdio.h>
float Camera::generate_ray(const CameraSample &sample, Ray *ray) const
{
    Point3f p_film(sample.p_film.x, sample.p_film.y, 0);
    Point3f p_camera = inverse(camera_to_screen) * raster_to_screen * p_film;

    // TODO: lerp time between shutter open/close
    *ray = Ray(Point3f(), normalize(Vector3f(p_camera)), INFINITY, sample.time);
    *ray = camera_to_world * (*ray);

    return 1;
}
