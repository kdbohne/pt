#include <iostream>
#include <fstream>
#include <sstream>
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

    inline T operator[](int i) const { assert((i >= 0) && (i <= 1)); return *(&x + i); }
    inline T &operator[](int i) { assert((i >= 0) && (i <= 1)); return *(&x + i); }

    T x, y;
};

template<typename T>
struct Vector3
{
    Vector3() : x(0), y(0), z(0) {}
    Vector3(T s) : x(s), y(s), z(s) {}
    Vector3(T x, T y, T z) : x(x), y(y), z(z) {}

    inline T operator[](int i) const { assert((i >= 0) && (i <= 2)); return *(&x + i); }
    inline T &operator[](int i) { assert((i >= 0) && (i <= 2)); return *(&x + i); }

    inline Vector3<T> operator+(const Vector3<T> &v) const { return Vector3<T>(x + v.x, y + v.y, z + v.z); }
    inline Vector3<T> operator-(const Vector3<T> &v) const { return Vector3<T>(x - v.x, y - v.y, z - v.z); }
    inline Vector3<T> &operator+=(const Vector3<T> &v) { x += v.x; y += v.y; z += v.z; return *this; }
    inline Vector3<T> &operator-=(const Vector3<T> &v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    inline Vector3<T> operator*(T s) const { return Vector3<T>(x * s, y * s, z * s); }
    inline Vector3<T> operator/(T s) const { return Vector3<T>(x / s, y / s, z / s); }
    inline Vector3<T> operator-() const { return Vector3<T>(-x, -y, -z); }

    inline float length_squared() const { return (x * x) + (y * y) + (z * z); }
    inline float length() const { return std::sqrt(length_squared()); }

    T x, y, z;
};

template<typename T>
inline Vector3<T> operator*(float s, const Vector3<T> &v) { return v * s; }

template<typename T>
inline Vector3<T> normalize(const Vector3<T> &v) { return v / v.length(); }

template<typename T>
inline Vector3<T> abs(const Vector3<T> &v)
{
    return Vector3<T>(std::abs(v.x), std::abs(v.y), std::abs(v.z));
}

template<typename T>
inline float dot(const Vector3<T> &a, const Vector3<T> &b)
{
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template<typename T>
inline float abs_dot(const Vector3<T> &a, const Vector3<T> &b)
{
    return std::abs(dot(a, b));
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

template<typename T>
inline float distance(const Vector3<T> &a, const Vector3<T> &b)
{
    return (a - b).length();
}

template<typename T>
inline float distance_squared(const Vector3<T> &a, const Vector3<T> &b)
{
    return (a - b).length_squared();
}

template<typename T>
inline int max_dimension(const Vector3<T> &v)
{
    return (v.x > v.y) ? ((v.x > v.z) ? 0 : 2) : ((v.y > v.z) ? 1 : 2);
}

template<typename T>
inline Vector3<T> permute(const Vector3<T> &v, int x, int y, int z)
{
    return Vector3<T>(v[x], v[y], v[z]);
}

template<typename T>
inline void coordinate_system(const Vector3<T> &v0, Vector3<T> *v1, Vector3<T> *v2)
{
    // TODO: does this work for right-handed coordinate systems?
    if (std::abs(v0.x) > std::abs(v0.y))
        *v1 = Vector3<T>(-v0.z, 0, v0.x) / std::sqrt(v0.x * v0.x + v0.z * v0.z);
    else
        *v1 = Vector3<T>(0, v0.z, -v0.y) / std::sqrt(v0.y * v0.y + v0.z * v0.z);

    *v2 = cross(v0, *v1);
}

typedef Vector2<float> Vector2f;
typedef Vector2<int>   Vector2i;
typedef Vector3<float> Vector3f;
typedef Vector3<int>   Vector3i;

inline Vector3f spherical_direction(float sin_theta, float cos_theta, float phi,
                                    const Vector3f &x, const Vector3f &y, const Vector3f &z)
{
    return sin_theta * std::cos(phi) * x +
           sin_theta * std::sin(phi) * y +
           cos_theta * z;
}

struct Ray
{
    Vector3f o;
    Vector3f d;
    mutable float tmax;

    float time;

    Ray() : o(Vector3f()), d(Vector3f()), tmax(INFINITY), time(0) {}
    Ray(const Vector3f &o, const Vector3f &d, float tmax, float time) : o(o), d(d), tmax(tmax), time(time) {}

    inline Vector3f evaluate(float t) const { return o + t * d; }
};

struct Matrix4x4
{
    float m[4][4];

    Matrix4x4() : Matrix4x4(1, 0, 0, 0,
                            0, 1, 0, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1)
    {
    }

    Matrix4x4(float m00, float m10, float m20, float m30,
              float m01, float m11, float m21, float m31,
              float m02, float m12, float m22, float m32,
              float m03, float m13, float m23, float m33)
    {
        m[0][0] = m00;
        m[0][1] = m01;
        m[0][2] = m02;
        m[0][3] = m03;

        m[1][0] = m10;
        m[1][1] = m11;
        m[1][2] = m12;
        m[1][3] = m13;

        m[2][0] = m20;
        m[2][1] = m21;
        m[2][2] = m22;
        m[2][3] = m23;

        m[3][0] = m30;
        m[3][1] = m31;
        m[3][2] = m32;
        m[3][3] = m33;
    }

    Matrix4x4(float mat[4][4]) { std::memcpy(m, mat, sizeof(float) * 16); }

    inline Matrix4x4 operator*(const Matrix4x4 &mat) const
    {
        Matrix4x4 r;
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                r.m[i][j] = m[0][j] * mat.m[i][0] +
                            m[1][j] * mat.m[i][1] +
                            m[2][j] * mat.m[i][2] +
                            m[3][j] * mat.m[i][3];
            }
        }
        return r;
    }
};

static Matrix4x4 transpose(const Matrix4x4 &m)
{
    return Matrix4x4(m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
                     m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
                     m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
                     m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
}

// NOTE: this is Inverse() from pbrt-v3.
static Matrix4x4 inverse(const Matrix4x4 &m)
{
    int indxc[4], indxr[4];
    int ipiv[4] = {0, 0, 0, 0};
    float minv[4][4];
    std::memcpy(minv, m.m, 4 * 4 * sizeof(float));
    for (int i = 0; i < 4; i++) {
        int irow = 0, icol = 0;
        float big = 0.f;
        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(minv[j][k]) >= big) {
                            big = float(std::abs(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1)
                        error("Singular matrix in MatrixInvert%s", "");
                }
            }
        }
        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (minv[icol][icol] == 0.f) error("Singular matrix in MatrixInvert%s", "");

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        float pivinv = 1. / minv[icol][icol];
        minv[icol][icol] = 1.;
        for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                float save = minv[j][icol];
                minv[j][icol] = 0;
                for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
            }
        }
    }
    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
        }
    }
    return Matrix4x4(minv);
}

// NOTE: this only exists to specify whether a Vector3f should be translated or
// not when being transformed by a Transform.
//     i.e. Vector3f v = ...;
//          Transform t = ...;
//          v = t * v;           <=>  t * vec4(v, 0)
//          v = t * Point3f(v);  <=>  t * vec4(v, 1)
template<typename T>
struct Point3
{
    Vector3<T> v;
    Point3<T>(const Vector3<T> &v) : v(v) {}
};

typedef Point3<float> Point3f;

struct Transform
{
    Matrix4x4 m;
    Matrix4x4 inv;

    Transform() {}
    Transform(const Matrix4x4 &m) : m(m), inv(inverse(m)) {}
    Transform(const Matrix4x4 &m, const Matrix4x4 &inv) : m(m), inv(inv) {}

    inline Transform operator*(const Transform &t) const
    {
        return Transform(m * t.m, t.inv * inv);
    }

    template<typename T>
    inline Vector3<T> operator*(const Vector3<T> &v) const
    {
        T x = v.x, y = v.y, z = v.z;
        return Vector3<T>(m.m[0][0] * x + m.m[1][0] * y + m.m[2][0] * z,
                          m.m[0][1] * x + m.m[1][1] * y + m.m[2][1] * z,
                          m.m[0][2] * x + m.m[1][2] * y + m.m[2][2] * z);
    }

    template<typename T>
    inline Vector3<T> operator*(const Point3<T> &p) const
    {
        // NOTE: assuming w = 1.
        T x = p.v.x, y = p.v.y, z = p.v.z;
        T xp = m.m[0][0] * x + m.m[1][0] * y + m.m[2][0] * z + m.m[3][0];
        T yp = m.m[0][1] * x + m.m[1][1] * y + m.m[2][1] * z + m.m[3][1];
        T zp = m.m[0][2] * x + m.m[1][2] * y + m.m[2][2] * z + m.m[3][2];
        T wp = m.m[0][3] * x + m.m[1][3] * y + m.m[2][3] * z + m.m[3][3];
        if (wp == 1)
            return Vector3<T>(xp, yp, zp);
        else
            return Vector3<T>(xp, yp, zp) / wp;
    }

    inline Ray operator*(const Ray &r) const
    {
        // TODO: handle error
        Vector3f o = (*this) * Point3f(r.o);
        Vector3f d = (*this) * r.d;

        return Ray(o, d, r.tmax, r.time);
    }
};

inline Transform inverse(const Transform &t)
{
    return Transform(t.inv, t.m);
}

static Transform translate(const Vector3f &delta)
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

static Transform scale(float x, float y, float z)
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

static Transform perspective(float vfov, float aspect, float n, float f)
{
    Matrix4x4 persp(1, 0,           0,                0,
                    0, 1,           0,                0,
                    0, 0, f / (f - n), -f * n / (f - n),
                    0, 0,           1,                0);

    float inv_tan = 1 / std::tan(vfov / 2);

    return scale(inv_tan / aspect, inv_tan, 1) * Transform(persp);
}

static Transform look_at(const Vector3f &eye, const Vector3f &target, const Vector3f &up)
{
    Vector3f forward = normalize(target - eye);
    Vector3f right = normalize(cross(forward, up));
    Vector3f local_up = normalize(cross(right, forward));

    Matrix4x4 camera_to_world;

    camera_to_world.m[0][0] = right.x;
    camera_to_world.m[0][1] = right.y;
    camera_to_world.m[0][2] = right.z;
    camera_to_world.m[0][3] = 0;

    camera_to_world.m[1][0] = local_up.x;
    camera_to_world.m[1][1] = local_up.y;
    camera_to_world.m[1][2] = local_up.z;
    camera_to_world.m[1][3] = 0;

    camera_to_world.m[2][0] = forward.x;
    camera_to_world.m[2][1] = forward.y;
    camera_to_world.m[2][2] = forward.z;
    camera_to_world.m[2][3] = 0;

    camera_to_world.m[3][0] = eye.x;
    camera_to_world.m[3][1] = eye.y;
    camera_to_world.m[3][2] = eye.z;
    camera_to_world.m[3][3] = 1;

    return camera_to_world;
}

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
    Transform camera_to_world;
    Transform camera_to_screen;

    Film *film;

    Camera(const Transform &camera_to_world, float vfov, Film *film)
        : camera_to_world(camera_to_world),
          camera_to_screen(perspective(vfov, (float)film->resolution.x / (float)film->resolution.y, 1e-2f, 1000.0f)),
          film(film)
    {
    }

    float generate_ray(const CameraSample &sample, Ray *ray) const
    {
        Vector3f screen_pos(sample.film_pos.x * 2 - 1, sample.film_pos.y * 2 - 1, 0);
        screen_pos = inverse(camera_to_screen) * screen_pos;
        screen_pos.z = 1;

        // TODO: lerp time between shutter open/close
        *ray = Ray(Vector3f(), normalize(screen_pos), INFINITY, sample.time);
        *ray = camera_to_world * (*ray);

        return 1;
    }
};

struct TriangleData
{
    int pi[3];
    int ni[3];
    int uvi[3];
};

struct MeshData
{
    std::vector<TriangleData> tris;
    std::vector<Vector3f> p;
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

                Vector3f p;
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
        for (Vector3f &p : this->data.p)
            p = transform * Point3f(p);
        for (Vector3f &n : this->data.n)
            n = transform * n;
    }
};

struct Entity;

struct Intersection
{
    Vector3f p;
    Vector3f n;
    float time;

    const Entity *entity;

    Ray spawn_ray(const Vector3f &d) const
    {
        return Ray(p, d, INFINITY, time);
    }
};

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

    virtual bool intersect(const Ray &ray, float *t, Intersection *intersection) const = 0;

    virtual Intersection sample(const Vector2f &u) const = 0;
    virtual Intersection sample(const Intersection &ref, const Vector2f &u) const = 0;

    virtual float area() const = 0;
    virtual float pdf(const Intersection &ref, const Vector3f &wi) const;
};

float Geometry::pdf(const Intersection &ref, const Vector3f &wi) const
{
    Ray ray = ref.spawn_ray(wi);

    float t;
    Intersection light_intersection;
    if (intersect(ray, &t, &light_intersection))
        return 0;

    float pdf = distance_squared(ref.p, light_intersection.p) / (abs_dot(light_intersection.n, -wi) * area());
    return pdf;
}

struct Light
{
    float emittance;

    Light(float emittance) : emittance(emittance) {}
};

struct Entity
{
    std::shared_ptr<Geometry> geometry;
    std::shared_ptr<Light> light;

    Entity(const std::shared_ptr<Geometry> &geometry, const std::shared_ptr<Light> &light) : geometry(geometry), light(light) {}

    bool intersect(const Ray &ray, float *t, Intersection *intersection) const
    {
        if (!geometry->intersect(ray, t, intersection))
            return false;

        intersection->entity = this;
        return true;
    }

    // TODO: return spectrum
    Vector3f radiance(const Intersection &intersection, const Vector3f &w) const
    {
        assert(light);
        return (dot(intersection.n, w) > 0) ? light->emittance : Vector3f(0);
    }

    // TODO: return spectrum
    // TODO: visibility tester
    Vector3f sample_light(const Intersection &ref, const Vector2f &u,
                          Vector3f *wi, float *pdf) const
    {
        assert(light);

        Intersection intersection = geometry->sample(ref, u);
        *wi = normalize(intersection.p - ref.p);
        *pdf = geometry->pdf(ref, *wi);
        // TODO: visibility tester
        return radiance(intersection, -(*wi));
    }
};

struct Sphere : public Geometry
{
    float radius;

    Sphere() : Sphere(Transform(), 0) {}
    Sphere(const Transform &object_to_world, float radius) : Geometry(object_to_world), radius(radius) {}

    float area() const override { return 4 * PI * (radius * radius); }

    bool intersect(const Ray &ray, float *t, Intersection *intersection) const override
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

        *t = (t0 > 0) ? t0 : t1;
        if (*t > r.tmax)
            return false;

        Vector3f hit = r.evaluate(*t);
        // TODO: refine hit point
        if ((hit.x == 0) && (hit.y == 0))
            hit.x = 1e-5f * radius;

        intersection->p = object_to_world * hit;
        intersection->n = normalize(intersection->p);

        ray.tmax = *t;

        return true;
    }

    Intersection sample(const Vector2f &u) const override
    {
        Vector3f obj = Vector3f(0, 0, 0) + radius * uniform_sample_sphere(u);

        Intersection i;
        i.n = normalize(object_to_world * obj);
        // TODO: reproject
        i.p = object_to_world * Point3f(obj);

        return i;
    }

    Intersection sample(const Intersection &ref, const Vector2f &u) const override
    {
        Vector3f center = object_to_world * Point3f(Vector3f(0, 0, 0));
        Vector3f wc = normalize(center - ref.p);
        Vector3f wcx, wcy;
        coordinate_system(wc, &wcx, &wcy);

        // TODO: offset_ray_origin()
        Vector3f origin = ref.p;
        if (distance_squared(origin, center) <= radius * radius)
            return sample(u);

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

        return i;
    }

    float pdf(const Intersection &ref, const Vector3f &wi) const override
    {
        Vector3f center = object_to_world * Point3f(Vector3f(0, 0, 0));
        // TODO: offset_ray_origin()
        Vector3f origin = ref.p;
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

    bool intersect(const Ray &ray, float *t, Intersection *intersection) const override
    {
        const Vector3f &p0 = mesh->data.p[pi[0]];
        const Vector3f &p1 = mesh->data.p[pi[1]];
        const Vector3f &p2 = mesh->data.p[pi[2]];

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
        *t = t_scaled * inv_det;

        Vector3f hit = b0 * p0 + b1 * p1 + b2 * p2;

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

        ray.tmax = *t;

        return true;
    }

    Intersection sample(const Vector2f &u) const override
    {
        // TODO FIXME
        return Intersection();
    }

    Intersection sample(const Intersection &ref, const Vector2f &u) const override
    {
        // TODO FIXME
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
            float t;
            if (entity.intersect(ray, &t, intersection))
                intersects = true;
        }

        return intersects;
    }

    void add_mesh(const std::shared_ptr<Mesh> &mesh)
    {
        for (const TriangleData &t : mesh->data.tris)
            entities.push_back(Entity(std::make_shared<Triangle>(mesh, t), nullptr));
    }
};

// TODO: return spectrum
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
    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(10, 10, 4)), 0.5), std::make_shared<Light>(1)));
    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(-1.25, 0, 0)), 0.1), std::make_shared<Light>(1)));
    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f(-3.75, 0, 0)), 0.03333), std::make_shared<Light>(1)));
    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 1.25, 0, 0)), 0.3), std::make_shared<Light>(1)));
    scene.entities.push_back(Entity(std::make_shared<Sphere>(translate(Vector3f( 3.75, 0, 0)), 0.9), std::make_shared<Light>(1)));

    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/floor.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate1.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate2.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate3.obj")));
    scene.add_mesh(std::make_shared<Mesh>(Transform(), load_obj("asset/veach_mi/plate4.obj")));

    Film film(Vector2i(768, 512));
    Camera camera(look_at(Vector3f(0, 2, 15), Vector3f(0, -2, 2.5), Vector3f(0, 1, 0)), radians(28), &film);

    render(scene, camera, film);
    film.write_ppm("out.ppm");

    return 0;
}
