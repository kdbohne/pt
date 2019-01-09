#include "camera.h"

Camera::Camera(const Transform &camera_to_world, float vfov, Film *film)
    : camera_to_world(camera_to_world),
      camera_to_screen(perspective(vfov, (float)film->resolution.x / (float)film->resolution.y, 1e-2, 1000.0)),
      film(film)
{
    float aspect_ratio = (float)film->resolution.x / (float)film->resolution.y;

    Bounds2f screen;
    if (aspect_ratio > 1)
    {
        screen.min.x = -aspect_ratio;
        screen.max.x =  aspect_ratio;
        screen.min.y = -1;
        screen.max.y =  1;
    }
    else
    {
        screen.min.x = -1;
        screen.max.x =  1;
        screen.min.y = -1 / aspect_ratio;
        screen.max.y =  1 / aspect_ratio;
    }

    Transform screen_to_raster;
    screen_to_raster = translate(Vector3f(-screen.min.x, -screen.min.y, 0));
    screen_to_raster = scale(1 / (screen.max.x - screen.min.x),
                             1 / (screen.max.y - screen.min.y),
                             1) * screen_to_raster;
    screen_to_raster = scale(film->resolution.x, film->resolution.y, 1) * screen_to_raster;

    raster_to_screen = inverse(screen_to_raster);
}

float Camera::generate_ray(const CameraSample &sample, Ray *ray) const
{
    Point3f p_film(sample.p_film.x, sample.p_film.y, 0);
    Point3f p_camera = inverse(camera_to_screen) * raster_to_screen * p_film;

    // TODO: lerp time between shutter open/close
    *ray = Ray(Point3f(), normalize(Vector3f(p_camera)), INFINITY, sample.time);
    *ray = camera_to_world * (*ray);

    return 1;
}
