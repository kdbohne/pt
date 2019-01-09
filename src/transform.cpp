#include "transform.h"
#include "math.h"

Transform Transform::operator*(const Transform &t) const
{
    return Transform(m * t.m, t.inv * inv);
}

Ray Transform::operator*(const Ray &r) const
{
    // TODO: handle error
    Point3f o = (*this) * Point3f(r.o);
    Vector3f d = (*this) * r.d;

    return Ray(o, d, r.tmax, r.time);
}

Intersection Transform::operator*(const Intersection &is) const
{
    Intersection result = is;
    // TODO: handle error
    result.p = (*this) * result.p;
    result.n = normalize((*this) * result.n);
    result.wo = normalize((*this) * result.wo);
    result.t = normalize((*this) * result.t);
    result.b = normalize((*this) * result.b);

    return result;
}

Transform inverse(const Transform &t)
{
    return Transform(t.inv, t.m);
}

Transform translate(const Vector3f &delta)
{
    Matrix4x4 m(1, 0, 0, delta.x,
                0, 1, 0, delta.y,
                0, 0, 1, delta.z,
                0, 0, 0, 1);

    Matrix4x4 inv(1, 0, 0, -delta.x,
                  0, 1, 0, -delta.y,
                  0, 0, 1, -delta.z,
                  0, 0, 0, 1);

    return Transform(m, inv);
}

Transform rotate(float angle, const Vector3f &axis)
{
    Vector3f a = normalize(axis);

    // TODO: should the angle parameter be specified in radians?
    float sin_theta = std::sin(radians(angle));
    float cos_theta = std::cos(radians(angle));

    Matrix4x4 m;
    m.m[0][0] = a.x * a.x + (1 - a.x * a.x) * cos_theta;
    m.m[0][1] = a.x * a.y * (1 - cos_theta) - a.z * sin_theta;
    m.m[0][2] = a.x * a.z * (1 - cos_theta) + a.y * sin_theta;

    m.m[1][0] = a.x * a.y * (1 - cos_theta) + a.z * sin_theta;
    m.m[1][1] = a.y * a.y + (1 - a.y * a.y) * cos_theta;
    m.m[1][2] = a.y * a.z * (1 - cos_theta) - a.x * sin_theta;

    m.m[2][0] = a.x * a.z * (1 - cos_theta) - a.y * sin_theta;
    m.m[2][1] = a.y * a.z * (1 - cos_theta) + a.x * sin_theta;
    m.m[2][2] = a.z * a.z + (1 - a.z * a.z) * cos_theta;

    return Transform(m, transpose(m));
}

Transform scale(float x, float y, float z)
{
    Matrix4x4 m(x, 0, 0, 0,
                0, y, 0, 0,
                0, 0, z, 0,
                0, 0, 0, 1);

    Matrix4x4 inv(1 / x,     0,     0, 0,
                      0, 1 / y,     0, 0,
                      0,     0, 1 / z, 0,
                      0,     0,     0, 1);

    return Transform(m, inv);
}

Transform perspective(float vfov, float aspect, float n, float f)
{
    Matrix4x4 persp(1, 0,           0,                0,
                    0, 1,           0,                0,
                    0, 0, f / (f - n), -f * n / (f - n),
                    0, 0,          -1,                0);

    float inv_tan = 1 / std::tan(vfov / 2);

    return scale(inv_tan, inv_tan, 1) * Transform(persp);
}

Transform look_at(const Point3f &eye, const Point3f &target, const Vector3f &up)
{
    Vector3f forward = normalize(target - eye);
    Vector3f right = normalize(cross(forward, up));
    Vector3f local_up = normalize(cross(right, forward)); // TODO: can this normalize() be removed?

    // TODO: check forward vs up to make sure they are not parallel

    Matrix4x4 camera_to_world;

    camera_to_world.m[0][0] = right.x;
    camera_to_world.m[0][1] = right.y;
    camera_to_world.m[0][2] = right.z;

    camera_to_world.m[1][0] = local_up.x;
    camera_to_world.m[1][1] = local_up.y;
    camera_to_world.m[1][2] = local_up.z;

    camera_to_world.m[2][0] = forward.x;
    camera_to_world.m[2][1] = forward.y;
    camera_to_world.m[2][2] = forward.z;

    camera_to_world.m[3][0] = eye.x;
    camera_to_world.m[3][1] = eye.y;
    camera_to_world.m[3][2] = eye.z;

    return Transform(inverse(camera_to_world), camera_to_world);
}
