#include "geometry.h"
#include "math.h"

#include <fstream>
#include <sstream>

Geometry::Geometry(const Transform &object_to_world, const Transform &world_to_object)
    : object_to_world(object_to_world), world_to_object(world_to_object)
{
}

Bounds3f Geometry::world_bounds() const
{
    return object_to_world * object_bounds();
}

Sphere::Sphere(const Transform &object_to_world, const Transform &world_to_object, float radius)
    : Geometry(object_to_world, world_to_object), radius(radius)
{
}

Bounds3f Sphere::object_bounds() const
{
    return Bounds3f(Point3f(-radius, -radius, -radius),
                    Point3f( radius,  radius,  radius));
}

bool Sphere::intersect(const Ray &ray, Intersection *intersection) const
{
    Ray r = inverse(object_to_world) * ray;
    float dx = r.d.x;
    float dy = r.d.y;
    float dz = r.d.z;
    float ox = r.o.x;
    float oy = r.o.y;
    float oz = r.o.z;

    float a = (dx * dx) + (dy * dy) + (dz * dz);
    float b = 2 * ((ox * dx) + (oy * dy) + (oz * dz));
    float c = (ox * ox) + (oy * oy) + (oz * oz) - (radius * radius);

    float t0, t1;
    if (!quadratic(a, b, c, &t0, &t1))
        return false;

    if ((t0 > r.tmax) || (t1 <= 0))
        return false;

    float t = (t0 > 0) ? t0 : t1;
    if (t > r.tmax)
        return false;

    Point3f hit = r.evaluate(t);
    // TODO: refine hit point
    if ((hit.x == 0) && (hit.y == 0))
        hit.x = 1e-5f * radius;

    Intersection object_is;
    object_is.p = hit;
    object_is.n = normalize(Normal3f(hit));
    object_is.wo = -r.d;
//    object_is.time = r.time;
//    object_is.uv = ; // TODO FIXME
//    object_is.dpdu = ; // TODO FIXME
//    object_is.dpdv = ; // TODO FIXME
//    object_is.t = ; // TODO FIXME
//    object_is.b = ; // TODO FIXME
//    object_is.entity = ; // TODO FIXME

    *intersection = object_to_world * object_is;

    ray.tmax = t;

    return true;
}

MeshData load_obj(const std::string &path)
{
    MeshData mesh;

    std::ifstream file(path);
    if (!file)
        error("Failed to load OBJ file: %s", path.c_str());

    std::string line;
    while (std::getline(file, line))
    {
        if (line[0] == 'v')
        {
            if (line[1] == ' ')
            {
                // Position
                std::istringstream ss(line.substr(2));

                Point3f p;
                ss >> p.x >> p.y >> p.z;
                mesh.p.push_back(p);
            }
            else if (line[1] == 'n')
            {
                std::istringstream ss(line.substr(3));

                Normal3f n;
                ss >> n.x >> n.y >> n.z;
                mesh.n.push_back(n);
            }
            else
            {
                // TODO: texture coordinates
                assert(false);
            }
        }
        else if (line[0] == 'f')
        {
            std::istringstream ss(line.substr(2));

            TriangleData t;
            for (int i = 0; i < 3; ++i)
            {
                int pi, ni = 0, uvi = 0;
                ss >> pi;
                if (ss.peek() == '/')
                {
                    ss.get();
                    if (ss.peek() != '/')
                        ss >> uvi;

                    ss.get();
                    ss >> ni;
                }

                // OBJ indexing starts at 1.
                t.pi[i] = pi - 1;
                t.ni[i] = ni - 1;
                t.uvi[i] = uvi - 1;
            }
            mesh.tris.push_back(t);
        }
    }

    return mesh;
}

Mesh::Mesh(const Transform &object_to_world, const Transform &world_to_object, const MeshData &data)
    : object_to_world(object_to_world), world_to_object(world_to_object), data(data)
{
    for (Point3f &p : this->data.p)
        p = object_to_world * p;
    for (Normal3f &n : this->data.n)
        n = object_to_world * n;
}

Triangle::Triangle(const Transform &object_to_world, const Transform &world_to_object, const Mesh *mesh, const TriangleData &data)
    : Geometry(object_to_world, world_to_object), mesh(mesh)
{
    pi[0] = data.pi[0];
    pi[1] = data.pi[1];
    pi[2] = data.pi[2];

    ni[0] = data.ni[0];
    ni[1] = data.ni[1];
    ni[2] = data.ni[2];

    uvi[0] = data.uvi[0];
    uvi[1] = data.uvi[1];
    uvi[2] = data.uvi[2];
}

Bounds3f Triangle::object_bounds() const
{
    const Point3f &p0 = mesh->data.p[pi[0]];
    const Point3f &p1 = mesh->data.p[pi[1]];
    const Point3f &p2 = mesh->data.p[pi[2]];

    return Union(Bounds3f(world_to_object * p0, world_to_object * p1), world_to_object * p2);
}

Bounds3f Triangle::world_bounds() const
{
    const Point3f &p0 = mesh->data.p[pi[0]];
    const Point3f &p1 = mesh->data.p[pi[1]];
    const Point3f &p2 = mesh->data.p[pi[2]];

    return Union(Bounds3f(p0, p1), p2);
}

bool Triangle::intersect(const Ray &ray, Intersection *intersection) const
{
    const Point3f &p0 = mesh->data.p[pi[0]];
    const Point3f &p1 = mesh->data.p[pi[1]];
    const Point3f &p2 = mesh->data.p[pi[2]];

    // Translate the triangle to the local ray coordinate system.
    Vector3f p0t = p0 - ray.o;
    Vector3f p1t = p1 - ray.o;
    Vector3f p2t = p2 - ray.o;

    // Determine which axis is largest, then permute the ray and
    // the triangle vertices so that this axis is now the z-axis.
    int kz = max_dimension(abs(ray.d));
    int kx = (kz + 1) % 3; // TODO: % vs if?
    int ky = (kx + 1) % 3; // TODO: % vs if?

    Vector3f d = permute(ray.d, kx, ky, kz);
    p0t = permute(p0t, kx, ky, kz);
    p1t = permute(p1t, kx, ky, kz);
    p2t = permute(p2t, kx, ky, kz);

    // Apply a shear transformation to the vertices to align them
    // with the ray facing in the +z direction.
    // NOTE: only applying x/y here; z is not needed for the next
    // two intersection tests.
    // TODO: precompute coefficients (sx, sy, sz), store in Ray
    float sx = -d.x / d.z;
    float sy = -d.y / d.z;
    float sz = 1 / d.z;
    p0t.x += sx * p0t.z;
    p0t.y += sy * p0t.z;
    p1t.x += sx * p1t.z;
    p1t.y += sy * p1t.z;
    p2t.x += sx * p2t.z;
    p2t.y += sy * p2t.z;

    // Signed edge function values.
    float e0 = p1t.x * p2t.y - p1t.y * p2t.x;
    float e1 = p2t.x * p0t.y - p2t.y * p0t.x;
    float e2 = p0t.x * p1t.y - p0t.y * p1t.x;

    // Fall back to double precision when testing triangle edges.
    if ((e0 == 0) || (e1 == 0) || (e2 == 0))
    {
        e0 = (float)((double)p1t.x * (double)p2t.y - (double)p1t.y * (double)p2t.x);
        e1 = (float)((double)p2t.x * (double)p0t.y - (double)p2t.y * (double)p0t.x);
        e2 = (float)((double)p0t.x * (double)p1t.y - (double)p0t.y * (double)p1t.x);
    }

    // Test triangle edges for containment.
    if (((e0 < 0) || (e1 < 0) || (e2 < 0)) && ((e0 > 0) || (e1 > 0) || (e2 > 0)))
        return false;

    // Ray is approaching edge-on.
    float e_det = e0 + e1 + e2;
    if (e_det == 0)
        return false;

    // Apply the z-component of the shear transformation now.
    p0t.z *= sz;
    p1t.z *= sz;
    p2t.z *= sz;

    // Scaled hit distance and corresponding range tests.
    float t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
    if ((e_det < 0) && ((t_scaled >= 0) || (t_scaled < ray.tmax * e_det)))
        return false;
    if ((e_det > 0) && ((t_scaled <= 0) || (t_scaled > ray.tmax * e_det)))
        return false;

    // Compute barycentric coordinates.
    float inv_e_det = 1 / e_det;
    float b0 = e0 * inv_e_det;
    float b1 = e1 * inv_e_det;
    float b2 = e2 * inv_e_det;
    float t = t_scaled * inv_e_det;

    Point2f uv[3];
    if (mesh->data.uv.empty())
    {
        uv[0] = Point2f(0, 0);
        uv[1] = Point2f(1, 0);
        uv[2] = Point2f(1, 1);
    }
    else
    {
        uv[0] = mesh->data.uv[uvi[0]];
        uv[1] = mesh->data.uv[uvi[1]];
        uv[2] = mesh->data.uv[uvi[2]];
    }

    // TODO: precompute, store in MeshData?
    Vector3f dpdu, dpdv;
    Vector2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
    Vector3f dp02 = p0 - p2, dp12 = p1 - p2;
    float uv_det = duv02[0] * duv12[1] - duv02[1] * duv12[0];
    if (uv_det == 0)
    {
        coordinate_system(normalize(cross(p2 - p0, p1 - p0)), &dpdu, &dpdv);
    }
    else
    {
        float inv_uv_det = 1 / uv_det;
        dpdu = ( duv12[1] * dp02 - duv02[1] * dp12) * inv_uv_det;
        dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inv_uv_det;
    }

    Point3f p_hit = b0 * p0 + b1 * p1 + b2 * p2;
    Point2f uv_hit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

    Normal3f n;
    if (mesh->data.n.empty() || (ni[0] == -1) || (ni[1] == -1) || (ni[2] == -1))
    {
        // TODO: precompute, store in MeshData?
        n = Normal3f(normalize(cross(dp02, dp12)));
    }
    else
    {
        const Normal3f &n0 = mesh->data.n[ni[0]];
        const Normal3f &n1 = mesh->data.n[ni[1]];
        const Normal3f &n2 = mesh->data.n[ni[2]];
        n = b0 * n0 + b1 * n1 + b2 * n2;
    }

    intersection->p = p_hit;
    intersection->n = n;
    intersection->wo = -ray.d;
//    intersection->uv = uv_hit;
//    intersection->dpdu = dpdu;
//    intersection->dpdv = dpdv;
//    intersection->dndu = dndu; // TODO
//    intersection->dndv = dndv; // TODO
    intersection->t = normalize(dpdu);
    intersection->b = cross(intersection->n, intersection->t);

    ray.tmax = t;

    return true;
}

#if 0
#include "geometry.h"
#include "math.h"
#include "sampling.h"

#include <fstream>
#include <sstream>

Intersection Geometry::sample(const Intersection &ref, const Point2f &u, float *pdf) const
{
    Intersection is = sample(u, pdf);

    Vector3f wi = is.p - ref.p;
    if (wi.length_squared() == 0)
    {
        *pdf = 0;
    }
    else
    {
        wi = normalize(wi);
        *pdf *= distance_squared(ref.p, is.p) / abs_dot(is.n, -wi);
    }

    return is;
}

float Geometry::pdf(const Intersection &ref, const Vector3f &wi) const
{
    Ray ray = ref.spawn_ray(wi);

    Intersection is;
    if (intersect(ray, &is))
        return 0;

    float pdf = distance_squared(ref.p, is.p) / (abs_dot(is.n, -wi) * area());
    return pdf;
}

Sphere::Sphere() : Sphere(Transform(), 0)
{
}

Sphere::Sphere(const Transform &object_to_world, float radius)
    : Geometry(object_to_world), radius(radius)
{
}

bool Sphere::intersect(const Ray &ray, Intersection *intersection) const
{
    Ray r = inverse(object_to_world) * ray;
    float dx = r.d.x;
    float dy = r.d.y;
    float dz = r.d.z;
    float ox = r.o.x;
    float oy = r.o.y;
    float oz = r.o.z;

    float a = (dx * dx) + (dy * dy) + (dz * dz);
    float b = 2 * ((ox * dx) + (oy * dy) + (oz * dz));
    float c = (ox * ox) + (oy * oy) + (oz * oz) - (radius * radius);

    float t0, t1;
    if (!quadratic(a, b, c, &t0, &t1))
        return false;

    if ((t0 > r.tmax) || (t1 <= 0))
        return false;

    float t = (t0 > 0) ? t0 : t1;
    if (t > r.tmax)
        return false;

    Point3f hit = r.evaluate(t);
    // TODO: refine hit point
    if ((hit.x == 0) && (hit.y == 0))
        hit.x = 1e-5f * radius;

    Intersection object_is;
    object_is.p = hit;
    object_is.n = normalize(Normal3f(hit));
    object_is.wo = -r.d;
    object_is.time = r.time;
//    object_is.uv = ; // TODO FIXME
//    object_is.dpdu = ; // TODO FIXME
//    object_is.dpdv = ; // TODO FIXME
//    object_is.t = ; // TODO FIXME
//    object_is.b = ; // TODO FIXME
//    object_is.entity = ; // TODO FIXME

    *intersection = object_to_world * object_is;

    ray.tmax = t;

    return true;
}

Intersection Sphere::sample(const Point2f &u, float *pdf) const
{
    Point3f obj = Point3f(0, 0, 0) + radius * uniform_sample_sphere(u);

    Intersection i;
    i.n = normalize(object_to_world * Normal3f(obj));

    // Reproject onto the surface of the sphere.
    obj *= radius / distance(obj, Point3f(0, 0, 0));
    i.p = object_to_world * obj;

    *pdf = 1 / area();

    return i;
}

Intersection Sphere::sample(const Intersection &ref, const Point2f &u, float *pdf) const
{
    Point3f center = object_to_world * Point3f(0, 0, 0);
    Vector3f wc = normalize(center - ref.p);
    Vector3f wcx, wcy;
    coordinate_system(wc, &wcx, &wcy);

    // TODO: offset_ray_origin()
    Point3f origin = ref.p;
    if (distance_squared(origin, center) <= radius * radius)
        return sample(u, pdf);

    float sin_theta_max2 = radius * radius / distance_squared(ref.p, center);
    float cos_theta_max = std::sqrt(std::max((float)0, 1 - sin_theta_max2));
    float cos_theta = (1 - u[0]) + u[0] * cos_theta_max;
    float sin_theta = std::sqrt(std::max((float)0, 1 - cos_theta * cos_theta));
    float phi = u[1] * 2 * PI;

    float dc = distance(ref.p, center);
    float ds = dc * cos_theta - std::sqrt(std::max((float)0, radius * radius - dc * dc * sin_theta * sin_theta));
    float cos_alpha = (dc * dc + radius * radius - ds * ds) / (2 * dc * radius);
    float sin_alpha = std::sqrt(std::max((float)0, 1 - cos_alpha * cos_alpha));

    Vector3f n_obj = spherical_direction(sin_alpha, cos_alpha, phi, -wcx, -wcy, -wc);
    Vector3f p_obj = radius * n_obj;

    Intersection i;
    // TODO: reproject
    i.p = object_to_world * Point3f(p_obj);
    i.n = object_to_world * Normal3f(n_obj);

    *pdf = uniform_cone_pdf(cos_theta_max);

    return i;
}

float Sphere::area() const
{
    return 4 * PI * (radius * radius);
}

float Sphere::pdf(const Intersection &ref, const Vector3f &wi) const
{
    Point3f center = object_to_world * Point3f(0, 0, 0);
    // TODO: offset_ray_origin()
    Point3f origin = ref.p;
    if (distance_squared(origin, center) <= radius * radius)
        return Geometry::pdf(ref, wi);

    float sin_theta_max2 = radius * radius / distance_squared(ref.p, center);
    float cos_theta_max = std::sqrt(std::max((float)0, 1 - sin_theta_max2));
    return uniform_cone_pdf(cos_theta_max);
}

Intersection Triangle::sample(const Point2f &u, float *pdf) const
{
    Point2f b = uniform_sample_triangle(u);

    const Point3f &p0 = mesh->data.p[pi[0]];
    const Point3f &p1 = mesh->data.p[pi[1]];
    const Point3f &p2 = mesh->data.p[pi[2]];

    Intersection is;
    is.p = b[0] * p0 + b[1] * p1 + (1 - b[0] - b[1]) * p2;
    if (mesh->data.n.empty())
    {
        // TODO: precompute, store in MeshData?
        is.n = Normal3f(normalize(cross(p1 - p0, p2 - p0)));
    }
    else
    {
        const Normal3f &n0 = mesh->data.n[ni[0]];
        const Normal3f &n1 = mesh->data.n[ni[1]];
        const Normal3f &n2 = mesh->data.n[ni[2]];
        is.n = b[0] * n0 + b[1] * n1 + (1 - b[0] - b[1]) * n2;
    }

    *pdf = 1 / area();

    return is;
}

float Triangle::area() const
{
    const Point3f &p0 = mesh->data.p[pi[0]];
    const Point3f &p1 = mesh->data.p[pi[1]];
    const Point3f &p2 = mesh->data.p[pi[2]];
    return 0.5 * cross(p1 - p0, p2 - p0).length();
}
#endif
