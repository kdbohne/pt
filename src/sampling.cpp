#include "sampling.h"
#include "math.h"

float uniform_float()
{
    // TODO: improve?
    return (float)drand48();
}

Vector3f uniform_sample_sphere(const Point2f &u)
{
    float z = 1 - 2 * u[0];
    float r = std::sqrt(std::max((float)0, (float)1 - z * z));
    float phi = 2 * PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

Point2f uniform_sample_triangle(const Point2f &u)
{
    float su0 = std::sqrt(u[0]);
    return Point2f(1 - su0, u[1] * su0);
}

Point2f concentric_sample_disk(const Point2f &u)
{
    Point2f u_offset = (float)2 * u - Vector2f(1, 1);
    if ((u_offset.x == 0) && (u_offset.y == 0))
        return Point2f(0, 0);

    float r, theta;
    if (std::abs(u_offset.x) > std::abs(u_offset.y))
    {
        r = u_offset.x;
        theta = PI_OVER_4 * (u_offset.y / u_offset.x);
    }
    else
    {
        r = u_offset.y;
        theta = PI_OVER_2 - PI_OVER_4 * (u_offset.x / u_offset.y);
    }

    return r * Point2f(std::cos(theta), std::sin(theta));
}

Vector3f cosine_sample_hemisphere(const Point2f &u)
{
    Point2f d = concentric_sample_disk(u);
    float z = std::sqrt(std::max((float)0, 1 - d.x * d.x - d.y * d.y));
    return Vector3f(d.x, d.y, z);
}

float uniform_cone_pdf(float cos_theta_max)
{
    return 1 / (2 * PI * (1 - cos_theta_max));
}
