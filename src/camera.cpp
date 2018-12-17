#include "camera.h"

float Camera::generate_ray(const CameraSample &sample, Ray *ray) const
{
    Point3f screen_pos(sample.film_pos.x * 2 - 1, sample.film_pos.y * 2 - 1, 0);
    screen_pos = inverse(camera_to_screen) * screen_pos;

    // TODO: lerp time between shutter open/close
    *ray = Ray(Point3f(), normalize(Vector3f(screen_pos)), INFINITY, sample.time);
    *ray = camera_to_world * (*ray);

    return 1;
}
