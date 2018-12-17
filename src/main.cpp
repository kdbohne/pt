#include "math.h"
#include "vector.h"
#include "point.h"
#include "ray.h"
#include "transform.h"
#include "camera.h"
#include "film.h"

#include <fstream>
#include <sstream>
#include <memory>

struct TriangleData
{
    int pi[3];
    int ni[3];
    int uvi[3];
};

struct MeshData
{
    std::vector<TriangleData> tris;
    std::vector<Point3f> p;
    std::vector<Vector3f> n;
};

static MeshData load_obj(const std::string &path)
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

                Vector3f n;
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

struct Mesh
{
    MeshData data;

    Mesh(const Transform &transform, const MeshData &data) : data(data)
    {
        for (Point3f &p : this->data.p)
            p = transform * Point3f(p);
        for (Vector3f &n : this->data.n)
            n = transform * n;
    }
};

struct Entity;

struct Intersection
{
    Point3f p;
    Vector3f n;
    Vector3f wo;

    float time;

    const Entity *entity;

    Intersection() {}
    Intersection(const Point3f &p, const Vector3f &n, const Vector3f &wo, float time)
        : p(p), n(n), wo(wo), time(time), entity(nullptr)
    {
    }

    Ray spawn_ray(const Vector3f &d) const
    {
        return Ray(p, d, INFINITY, time);
    }
};

// TODO: make member of Transform
inline Intersection operator*(const Transform &t, const Intersection &in)
{
    Intersection result;
    // TODO: handle error
    result.p = t * Point3f(in.p);
    result.n = t * in.n;
    result.wo = t * in.wo;
    result.time = in.time;
    result.entity = in.entity;

    return result;
}

static Vector3f uniform_sample_sphere(const Vector2f &u)
{
    float z = 1 - 2 * u[0];
    float r = std::sqrt(std::max((float)0, (float)1 - z * z));
    float phi = 2 * PI * u[1];
    return Vector3f(r * std::cos(phi), r * std::sin(phi), z);
}

static float uniform_cone_pdf(float cos_theta_max)
{
    return 1 / (2 * PI * (1 - cos_theta_max));
}

struct Geometry
{
    Transform object_to_world; // TODO: cache?

    Geometry(const Transform &object_to_world) : object_to_world(object_to_world) {}
    virtual ~Geometry() {}

    virtual bool intersect(const Ray &ray, Intersection *intersection) const = 0;

    virtual Intersection sample(const Vector2f &u, float *pdf) const = 0;
    virtual Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const = 0;

    virtual float area() const = 0;
    virtual float pdf(const Intersection &ref, const Vector3f &wi) const;
};

float Geometry::pdf(const Intersection &ref, const Vector3f &wi) const
{
    Ray ray = ref.spawn_ray(wi);

    Intersection light_intersection;
    if (intersect(ray, &light_intersection))
        return 0;

    float pdf = distance_squared(ref.p, light_intersection.p) / (abs_dot(light_intersection.n, -wi) * area());
    return pdf;
}

struct Light
{
    float emittance;
    std::shared_ptr<Geometry> geometry;

    Light(float emittance, const std::shared_ptr<Geometry> &geometry) : emittance(emittance), geometry(geometry) {}

    // TODO: return spectrum
    Vector3f radiance(const Intersection &intersection, const Vector3f &w) const
    {
        return (dot(intersection.n, w) > 0) ? emittance : Vector3f(0);
    }

    // TODO: return spectrum
    // TODO: visibility tester
    Vector3f sample_incident_radiance(const Intersection &ref, const Vector2f &u,
                                      Vector3f *wi, float *pdf) const
    {
        Intersection geo_inter = geometry->sample(ref, u, pdf);
        *wi = normalize(geo_inter.p - ref.p);
        // TODO: visibility tester
        return radiance(geo_inter, -(*wi));
    }
};

struct Entity
{
    std::shared_ptr<Geometry> geometry;
    std::shared_ptr<Light> light;

    Entity(const std::shared_ptr<Geometry> &geometry, const std::shared_ptr<Light> &light) : geometry(geometry), light(light) {}

    bool intersect(const Ray &ray, Intersection *intersection) const
    {
        if (!geometry->intersect(ray, intersection))
            return false;

        intersection->entity = this;
        return true;
    }

#if 0
    // TODO: return spectrum
    // TODO: visibility tester
    Vector3f sample(const Intersection &ref, const Vector2f &u,
                          Vector3f *wi, float *pdf) const
    {
        assert(light);

        Intersection intersection = geometry->sample(ref, u);
        *wi = normalize(intersection.p - ref.p);
        *pdf = geometry->pdf(ref, *wi);
        // TODO: visibility tester
        return light->radiance(intersection, -(*wi));
    }
#endif
};

struct Sphere : public Geometry
{
    float radius;

    Sphere() : Sphere(Transform(), 0) {}
    Sphere(const Transform &object_to_world, float radius) : Geometry(object_to_world), radius(radius) {}

    float area() const override { return 4 * PI * (radius * radius); }

    bool intersect(const Ray &ray, Intersection *intersection) const override
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

        *intersection = object_to_world * Intersection(hit, normalize(Vector3f(hit)), -r.d, r.time);

        ray.tmax = t;

        return true;
    }

    Intersection sample(const Vector2f &u, float *pdf) const override
    {
        Point3f obj = Point3f(0, 0, 0) + radius * uniform_sample_sphere(u);

        Intersection i;
        i.n = normalize(object_to_world * Vector3f(obj));
        // TODO: reproject
        i.p = object_to_world * obj;

        *pdf = 1 / area();

        return i;
    }

    Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const override
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
        i.n = object_to_world * n_obj;

        *pdf = uniform_cone_pdf(cos_theta_max);

        return i;
    }

    float pdf(const Intersection &ref, const Vector3f &wi) const override
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
};

struct Triangle : public Geometry
{
    Triangle(const std::shared_ptr<Mesh> &mesh, const TriangleData &data) : Geometry(Transform()), mesh(mesh)
    {
        pi[0] = data.pi[0];
        pi[1] = data.pi[1];
        pi[2] = data.pi[2];

        // TODO FIXME
        ni[0] = data.ni[0];
        ni[1] = data.ni[1];
        ni[2] = data.ni[2];
    }

    std::shared_ptr<Mesh> mesh;
    // TODO: reduce memory usage here
    int pi[3];
    int ni[3];

    bool intersect(const Ray &ray, Intersection *intersection) const override
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
        float sz = 1.0f / d.z;
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
        float det = e0 + e1 + e2;
        if (det == 0)
            return false;

        // Apply the z-component of the shear transformation now.
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;

        // Scaled hit distance and corresponding range tests.
        float t_scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if ((det < 0) && ((t_scaled >= 0) || (t_scaled < ray.tmax * det)))
            return false;
        if ((det > 0) && ((t_scaled <= 0) || (t_scaled > ray.tmax * det)))
            return false;

        // Compute barycentric coordinates.
        float inv_det = 1 / det;
        float b0 = e0 * inv_det;
        float b1 = e1 * inv_det;
        float b2 = e2 * inv_det;
        float t = t_scaled * inv_det;

        Point3f hit = b0 * p0 + b1 * p1 + b2 * p2;

        Vector3f n;
        if ((ni[0] == -1) || (ni[1] == -1) || (ni[2] == -1))
        {
            // TODO: precompute, store in MeshData?
            Vector3f dp02 = p0 - p2;
            Vector3f dp12 = p1 - p2;
            n = normalize(cross(dp02, dp12));
        }
        else
        {
            const Vector3f &n0 = mesh->data.n[ni[0]];
            const Vector3f &n1 = mesh->data.n[ni[1]];
            const Vector3f &n2 = mesh->data.n[ni[2]];
            n = b0 * n0 + b1 * n1 + b2 * n2;
        }

        intersection->p = hit;
        intersection->n = n;

        ray.tmax = t;

        return true;
    }

    Intersection sample(const Vector2f &u, float *pdf) const override
    {
        // TODO FIXME
        *pdf = 0;
        return Intersection();
    }

    Intersection sample(const Intersection &ref, const Vector2f &u, float *pdf) const override
    {
        // TODO FIXME
        *pdf = 0;
        return Intersection();
    }

    float area() const override
    {
        // TODO FIXME
        return 0;
    }

    float pdf(const Intersection &ref, const Vector3f &wi) const override
    {
        // TODO FIXME
        return 0;
    }
};

struct Scene
{
    std::vector<Entity> entities;

    bool intersect(const Ray &ray, Intersection *intersection) const
    {
        bool intersects = false;

        for (const Entity &entity : entities)
        {
            // TODO CLEANUP: get rid of t out parameter?
            if (entity.intersect(ray, intersection))
                intersects = true;
        }

        return intersects;
    }

    void add_light(const std::shared_ptr<Geometry> geo, float emittance)
    {
        entities.push_back(Entity(geo, std::make_shared<Light>(emittance, geo)));
    }

    void add_mesh(const std::shared_ptr<Mesh> &mesh)
    {
        for (const TriangleData &t : mesh->data.tris)
            entities.push_back(Entity(std::make_shared<Triangle>(mesh, t), nullptr));
    }
};

static float uniform_float()
{
    // TODO: improve
    return (float)drand48();
}

static float power_heuristic(int nf, float fpdf, int ng, float gpdf)
{
    float f = nf * fpdf;
    float g = ng * gpdf;
    return (f * f) / (f * f + g * g);
}

// TODO: spectrum
static Vector3f estimate_direct(const Light &light, const Intersection &intersection)
{
    // TODO: spectrum
    Vector3f direct(0);

    Vector2f u_light(uniform_float(), uniform_float());
    Vector3f wi;
    float light_pdf = 0;
    float scattering_pdf = 0;
    Vector3f incident = light.sample_incident_radiance(intersection, u_light, &wi, &light_pdf);
    if ((light_pdf > 0)/* && !incident.is_black()*/) // TODO: spectrum
    {
        Vector3f f(1); // TODO FIXME
        if (true/*!f.is_black()*/) // TODO FIXME
        {
            // TODO: scattering PDF
            // TODO: visibility
            float weight = power_heuristic(1, light_pdf, 1, scattering_pdf);
            direct += Vector3f(f.x * incident.x, f.y * incident.y, f.z * incident.z) * weight / light_pdf;
//            direct += f * incident * weight / light_pdf;
        }
    }

    return direct;
}

// TODO: spectrum
static Vector3f uniform_sample_all_lights(const Scene &scene, const Intersection &intersection)
{
    Vector3f radiance(0);

    // TODO: keep separate list of lights?
    for (const Entity &entity : scene.entities)
    {
        const std::shared_ptr<Light> &light = entity.light;
        if (!light)
            continue;

        const int samples_count = 1;

        Vector3f direct(0);
        for (int i = 0; i < samples_count; ++i)
            direct += estimate_direct(*light, intersection);

        radiance += direct / (float)samples_count;
    }

    return radiance;
}

// TODO: return spectrum
static Vector3f incident_radiance(const Scene &scene, const Ray &ray)
{
    Intersection intersection;
    if (!scene.intersect(ray, &intersection))
    {
        // TODO: light emission
        return Vector3f(0, 0, 0);
    }

    Vector3f radiance;

    const std::shared_ptr<Light> &light = intersection.entity->light;
    if (light)
        radiance += light->radiance(intersection, intersection.wo);

    radiance += uniform_sample_all_lights(scene, intersection);

    return radiance;
//    return intersection.n * 0.5 + Vector3f(0.5, 0.5, 0.5);
}

static void render(const Scene &scene, const Camera &camera, Film &film)
{
    for (int y = 0; y < film.resolution.y; ++y)
    {
        for (int x = 0; x < film.resolution.x; ++x)
        {
            CameraSample cs;
            cs.film_pos.x = (float)x / (float)(film.resolution.x - 1);
            cs.film_pos.y = 1 - (float)y / (float)(film.resolution.y - 1);
            cs.time = 0; // TODO

            Ray ray;
            camera.generate_ray(cs, &ray);

            Vector3f radiance = incident_radiance(scene, ray);
            film.pixels[y * film.resolution.x + x].rgb = radiance;
        }
    }
}

int main(int argc, char *argv[])
{
    Scene scene;
#if 1
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(10, 10, 4)), 0.5), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-1.25, 0, 0)), 0.1), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f(-3.75, 0, 0)), 0.03333), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 1.25, 0, 0)), 0.3), 1);
    scene.add_light(std::make_shared<Sphere>(translate(Vector3f( 3.75, 0, 0)), 0.9), 1);

    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/floor.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate1.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate2.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate3.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate4.obj")));

    Film film(Vector2i(768, 512) / 2);
    Camera camera(look_at(Vector3f(0, 2, 15), Vector3f(0, -2, 2.5), Vector3f(0, 1, 0)), radians(28), &film);
#endif

#if 0
    Film film(Vector2i(640, 360));
    Camera camera(look_at(Vector3f(8, 2, 3), Vector3f(0, 0, 0), Vector3f(0, 1, 0)), radians(40), &film);

    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(0, 0, 0)), 1), nullptr));
    {
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(0, -1000, 0)), 1000), nullptr));

        for (int a = -11; a < 11; ++a)
        {
            for (int b = -11; b < 11; ++b)
            {
                float choice = drand48();
                Vector3f center(a + 0.9 * drand48(), 0.2, b + 0.9 * drand48());
                if ((center - Vector3f(4.0, 0.2, 0.0)).length() > 0.9)
                {
                    if (choice < 0.8)
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                    else if (choice < 0.95)
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                    else
                        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(center), 0.2), nullptr));
                }
            }
        }

        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 0, 1, 0)), 1), nullptr));
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(-4, 1, 0)), 1), nullptr));
        scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 4, 1, 0)), 1), nullptr));
    }
#endif

    render(scene, camera, film);
    film.write_ppm("out.ppm");

    return 0;
}
