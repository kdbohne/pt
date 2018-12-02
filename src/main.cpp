#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>
#include <limits>
#include <cmath>
#include <cstring>

#define UNUSED(x) ((void)(x))

void report_(const char *prefix, bool fatal, const char *format, ...) __attribute__((__format__(__printf__, 3, 4)));
void report_(const char *prefix, bool fatal, const char *format, ...)
{
    std::printf("[%s] ", prefix);

    va_list args;
    va_start(args, format);
    std::vprintf(format, args); // TODO LOG
    va_end(args);

    std::printf("\n");

    if (fatal)
        __builtin_trap();
}

#define error(format, ...) report_("ERROR", false, format, __VA_ARGS__)
#define fatal(format, ...) report_("FATAL", true, format, __VA_ARGS__)
#define assert(expr) \
    ((expr) ? (void)0 : fatal("(%s:%d) Assertion failed: %s", __FILE__, __LINE__, #expr))

#undef INFINITY
static constexpr float INFINITY = std::numeric_limits<float>::infinity();
static constexpr float PI = 3.14159265358979323846;

inline float radians(float deg)
{
    return deg * (PI / 180.0f);
}

// Numerically-stable quadratic equation as written in Physically Based
// Rendering, Third Edition, pp. 1079-1080.
inline bool quadratic(float a, float b, float c, float *t0, float *t1)
{
    double discrim = (double)b * (double)b - 4 * (double)a * (double)c;
    if (discrim < 0)
        return false;

    double root_discrim = std::sqrt(discrim);

    double q;
    if (b < 0)
        q = -0.5 * (b - root_discrim);
    else
        q = -0.5 * (b + root_discrim);

    *t0 = q / a;
    *t1 = c / q;

    if (*t0 > *t1)
        std::swap(*t0, *t1);

    return true;
}

template<typename T>
struct Vector2
{
    Vector2() : x(0), y(0) {}
    Vector2(T x, T y) : x(x), y(y) {}

    T x, y;
};

template<typename T>
struct Vector3
{
    Vector3() : x(0), y(0), z(0) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    inline T operator[](int i) const { return *(&x + i); }
    inline T &operator[](int i) { return *(&x + i); }

    inline Vector3<T> operator+(const Vector3<T> &v) const { return Vector3<T>(x + v.x, y + v.y, z + v.z); }
    inline Vector3<T> operator-(const Vector3<T> &v) const { return Vector3<T>(x - v.x, y - v.y, z - v.z); }
    inline Vector3<T> &operator+=(const Vector3<T> &v) { x += v.x; y += v.y; z += v.z; return *this; }
    inline Vector3<T> &operator-=(const Vector3<T> &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    inline Vector3<T> operator*(T s) const { return Vector3<T>(x * s, y * s, z * s); }
    inline Vector3<T> operator/(T s) const { return Vector3<T>(x / s, y / s, z / s); }

    inline float length_squared() const { return (x * x) + (y * y) + (z * z); }
    inline float length() const { return std::sqrt(length_squared()); }

    T x, y, z;
};

template<typename T>
inline Vector3<T> operator*(float s, const Vector3<T> &v) { return v * s; }

template<typename T>
inline Vector3<T> normalize(const Vector3<T> &v) { return v / v.length(); }

template<typename T>
inline float dot(const Vector3<T> &a, const Vector3<T> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template<typename T>
inline Vector3<T> cross(const Vector3<T> &a, const Vector3<T> &b)
{
    double ax = a.x, ay = a.y, az = a.z;
    double bx = b.x, by = b.y, bz = b.z;
    return Vector3<T>((ay * bz) - (az * by),
                      (az * bx) - (ax * bz),
                      (ax * by) - (ay * bx));
}

typedef Vector2<float> Vector2f;
typedef Vector2<int>   Vector2i;
typedef Vector3<float> Vector3f;
typedef Vector3<int>   Vector3i;

struct Ray
{
    Vector3f o;
    Vector3f d;
    mutable float tmax;

    float time;

    inline Vector3f evaluate(float t) const { return o + t * d; }
};

struct Pixel
{
    // TODO: spectrum
    Vector3f rgb;
};

struct Film
{
    Vector2i resolution;
    std::vector<Pixel> pixels;

    Film(Vector2i resolution) : resolution(resolution), pixels(resolution.x * resolution.y) {}

    void write_ppm(const std::string &path) const
    {
        std::ofstream file(path);

        file << "P3\n";
        file << resolution.x << " " << resolution.y << "\n";
        file << "255" << "\n";

        for (int y = 0; y < resolution.y; ++y)
        {
            for (int x = 0; x < resolution.x; ++x)
            {
                // TODO: tone mapping
                const Pixel *p = &pixels[y * resolution.x + x];
                int r = (int)(p->rgb[0] * 255);
                int g = (int)(p->rgb[1] * 255);
                int b = (int)(p->rgb[2] * 255);

                file << r << " " << g << " " << b << " ";
            }

            file << "\n";
        }
    }
};

struct CameraSample
{
    Vector2f film_pos;
    float time;
};

struct Camera
{
    Vector3f pos;
    float vfov;

    Film *film;

    float half_height;
    float half_width;

    Camera(Vector3f pos, float vfov, Film *film) : pos(pos), vfov(vfov)
    {
        float aspect_ratio = (float)film->resolution.x / (float)film->resolution.y;
        half_height = std::tan(vfov / 2);
        half_width = half_height * aspect_ratio;
    }

    float generate_ray(const CameraSample &sample, Ray *ray) const
    {
        Vector3f forward(0, 0, -1);
        Vector3f right = normalize(cross(forward, Vector3f(0, 1, 0)));
        Vector3f up = normalize(cross(right, forward));

        float nx = 2 * sample.film_pos.x - 1;
        float ny = 2 * sample.film_pos.y - 1;
        ray->o = pos;
        ray->d = normalize(forward + (right * nx * half_width) + (up * ny * half_height));

        ray->tmax = INFINITY;
        ray->time = sample.time;

        return 1;
    }
};

struct Intersection
{
    Vector3f p;
    Vector3f n;
};

struct Sphere
{
    Sphere() : Sphere(Vector3f(), 0) {}
    Sphere(Vector3f center, float radius) : center(center), radius(radius) {}

    Vector3f center;
    float radius;
};

struct Scene
{
    std::vector<Sphere> spheres;

    bool intersect(const Ray &ray, Intersection *intersection) const
    {
        bool intersects = false;
        for (const Sphere &sphere : spheres)
        {
            // TODO: replace with transform
            Ray r = ray;
            r.o += sphere.center;

            // x^2 + y^2 + z^2 - r^2 = 0
            // (ox + t*dx)^2 + (oy + t*dy)^2 + (oz + t*dz)^2 - r^2 = 0
            // ox^2 + 2*ox*t*dx + t^2*dx^2 + oy^2 + 2*oy*t*dy + t^2*dy^2 + oz^2 + 2*oz*t*dz + t^2*dz^2 - r^2 = 0
            // ox^2 + oy^2 + oz^2 - r^2 + 2*ox*t*dx + 2*oy*t*dy + 2*oz*t*dz + t^2*dx^2 + t^2*dy^2 + t^2*dz^2 = 0
            // (ox^2 + oy^2 + oz^2 - r^2) + t*2*(ox*dx + oy*dy + oz*dz) + t^2*(dx^2 + dy^2 + dz^2) = 0
            //
            // a*t^2 + b*t + c = 0
            // a = dx^2 + dy^2 + dz^2
            // b = 2*(ox*dx + oy*dy + oz*dz)
            // c = dx^2 + dy^2 + dz^2
            float dx = r.d.x;
            float dy = r.d.y;
            float dz = r.d.z;
            float ox = r.o.x;
            float oy = r.o.y;
            float oz = r.o.z;

            float a = (dx * dx) + (dy * dy) + (dz * dz);
            float b = 2 * ((ox * dx) + (oy * dy) + (oz * dz));
            float c = (ox * ox) + (oy * oy) + (oz * oz) - (sphere.radius * sphere.radius);

            float t0, t1;
            if (!quadratic(a, b, c, &t0, &t1))
                continue;

            if ((t0 > r.tmax) || (t1 <= 0))
                continue;

            float t = (t0 > 0) ? t0 : t1;
            if (t > r.tmax)
                continue;

            ray.tmax = t;

            Vector3f hit = r.evaluate(t);
            // TODO: refine hit point
            if ((hit.x == 0) && (hit.y == 0))
                hit.x = 1e-5f * sphere.radius;

            intersection->p = hit - sphere.center;
            intersection->n = normalize(intersection->p);
            
            // TODO: just do return (intersection->p != Vector3f()) or something?
            intersects = true;
        }

        return intersects;
    }
};

static Vector3f incident_radiance(const Scene &scene, const Ray &ray)
{
    Intersection intersection;
    if (!scene.intersect(ray, &intersection))
    {
        // TODO: light emission
        return Vector3f(0, 0, 0);
    }

    return intersection.n * 0.5 + Vector3f(0.5, 0.5, 0.5);
}

static void render(const Scene &scene, const Camera &camera, Film &film)
{
    for (int y = 0; y < film.resolution.y; ++y)
    {
        for (int x = 0; x < film.resolution.x; ++x)
        {
            CameraSample cs;
            cs.film_pos.x = (float)x / (float)(film.resolution.x - 1);
            cs.film_pos.y = (float)y / (float)(film.resolution.y - 1);
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
    UNUSED(argc);
    UNUSED(argv);

    Scene scene;
    scene.spheres.push_back(Sphere(Vector3f(), 1));
    scene.spheres.push_back(Sphere(Vector3f(0, -1000, 0), 1000));

    Film film(Vector2i(300, 200));
    Camera camera(Vector3f(0, 0, 5), radians(60), &film);

    render(scene, camera, film);
    film.write_ppm("out.ppm");

    return 0;
}
