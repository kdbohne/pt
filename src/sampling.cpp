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

Distribution1d::Distribution1d(const float *f, int n)
    : func(f, f + n), cdf(n + 1)
{
    cdf[0] = 0;
    for (int i = 1; i < n + 1; ++i)
        cdf[i] = cdf[i - 1] + func[i - 1] / n;

    func_int = cdf[n];
    if (func_int == 0)
    {
        for (int i = 1; i < n + 1; ++i)
            cdf[i] = (float)i / (float)n;
    }
    else
    {
        for (int i = 1; i < n + 1; ++i)
            cdf[i] /= func_int;
    }
}

float Distribution1d::sample_continuous(float u, float *pdf, int *offset) const
{
    int off = find_interval((int)cdf.size(), [&](int index) { return cdf[index] <= u; });
    if (offset)
        *offset = off;

    float du = u - cdf[off];
    if ((cdf[off + 1] - cdf[off]) > 0)
    {
        assert(cdf[off + 1] > cdf[off]);
        du /= (cdf[off + 1] - cdf[off]);
    }
    assert(!std::isnan(du));

    if (pdf)
        *pdf = (func_int > 0) ? func[off] / func_int : 0;

    int count = (int)func.size();
    return (off + du) / count;
}

Distribution2d::Distribution2d(float *func, int nu, int nv)
{
    p_conditional_v.reserve(nv);
    for (int v = 0; v < nv; ++v)
        p_conditional_v.emplace_back(new Distribution1d(&func[v * nu], nu));

    std::vector<float> marginal_func;
    marginal_func.reserve(nv);
    for (int v = 0; v < nv; ++v)
        marginal_func.push_back(p_conditional_v[v]->func_int);

    p_marginal = new Distribution1d(&marginal_func[0], nv);
}

Point2f Distribution2d::sample_continuous(const Point2f &u, float *pdf) const
{
    float pdfs[2];
    int v;
    float d1 = p_marginal->sample_continuous(u[1], &pdfs[1], &v);
    float d0 = p_conditional_v[v]->sample_continuous(u[0], &pdfs[0]);

    *pdf = pdfs[0] * pdfs[1];

    return Point2f(d0, d1);
}
