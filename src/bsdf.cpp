#include "bsdf.h"
#include "math.h"
#include "sampling.h"
#include "intersection.h"

bool Bxdf::matches(BxdfType t) const
{
    return (type & t) == type;
}

Spectrum Bxdf::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf) const
{
    *wi = cosine_sample_hemisphere(u);
    if (wo.z < 0)
        wi->z *= -1;

    *pdf = this->pdf(wo, *wi);

    return f(wo, *wi);
}

float Bxdf::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    if (!same_hemisphere(wo, wi))
        return 0;

    return abs_cos_theta(wi) * INV_PI;
}

void Bsdf::add(Bxdf *bxdf)
{
    assert(bxdfs_count < MAX_BXDFS_COUNT);
    bxdfs[bxdfs_count++] = bxdf;
}

int Bsdf::get_matching_components_count(BxdfType flags) const
{
    int count = 0;

    for (int i = 0; i < bxdfs_count; ++i)
    {
        if (bxdfs[i]->matches(flags))
            ++count;
    }

    return count;
}

static Vector3f world_to_local(const Intersection &its, const Vector3f &v)
{
    return Vector3f(dot(v, its.t), dot(v, its.b), dot(v, its.n));
}

static Vector3f local_to_world(const Intersection &its, const Vector3f &v)
{
    return Vector3f(its.t.x * v.x + its.b.x * v.y + its.n.x * v.z,
                    its.t.y * v.x + its.b.y * v.y + its.n.y * v.z,
                    its.t.z * v.x + its.b.z * v.y + its.n.z * v.z);
}

Spectrum Bsdf::f(const Intersection &its, const Vector3f &wo, const Vector3f &wi, BxdfType flags) const
{
    Vector3f wo_local = world_to_local(its, wo);
    if (wo_local.z == 0)
        return Spectrum(0);
    Vector3f wi_local = world_to_local(its, wi);

    bool reflect = (dot(wi, its.n) * dot(wo, its.n) > 0);

    Spectrum f(0);
    for (int i = 0; i < bxdfs_count; ++i)
    {
        Bxdf *bxdf = bxdfs[i];
        if (bxdf->matches(flags))
        {
            if ((reflect && (bxdf->type & BSDF_REFLECTION)) ||
                (!reflect && (bxdf->type & BSDF_TRANSMISSION)))
            {
                    f += bxdf->f(wo_local, wi_local);
            }
        }
    }

    return f;
}

Spectrum Bsdf::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType type, BxdfType *sampled_type) const
{
    int components_count = get_matching_components_count(type);
    if (components_count == 0)
    {
        *pdf = 0;
        *sampled_type = (BxdfType)0;
        return Spectrum(0);
    }

    int component = std::min((int)std::floor(u[0] * components_count), components_count - 1);

    // TODO CLEANUP: ugly iteration
    Bxdf *bxdf = nullptr;
    int count = components_count;
    for (int i = 0; i < bxdfs_count; ++i)
    {
        if (bxdfs[i]->matches(type) && (count-- == 0))
        {
            bxdf = bxdfs[i];
            break;
        }
    }
    assert(bxdf);

    Point2f u_remapped(u[0] * components_count - component, u[1]);

    if (sampled_type)
        *sampled_type = bxdf->type;

    Vector3f wo_local = world_to_local(its, wo);
    Vector3f wi_local;
    *pdf = 0;
    Spectrum f = bxdf->sample_f(its, wo_local, &wi_local, u_remapped, pdf);
    if (*pdf == 0)
        return 0;

    *wi = local_to_world(its, wi_local);

    if (!(bxdf->type & BSDF_SPECULAR) && (components_count > 1))
    {
        for (int i = 0; i < bxdfs_count; ++i)
        {
            if ((bxdfs[i] != bxdf) && (bxdfs[i]->matches(type)))
                *pdf += bxdfs[i]->pdf(wo_local, wi_local);
        }
    }

    if (components_count > 1)
        *pdf /= components_count;

    if (!(bxdf->type & BSDF_SPECULAR) && (components_count > 1))
    {
        bool reflect = (dot(*wi, its.n) * dot(wo, its.n) > 0);

        f = 0;
        for (int i = 0; i < bxdfs_count; ++i)
        {
            if (bxdfs[i]->matches(type))
            {
                if (reflect && (bxdfs[i]->type & BSDF_REFLECTION))
                {
                    if (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))
                        f += bxdfs[i]->f(wo_local, wi_local);
                }
            }
        }
    }

    return f;
}

#if 0
float Bsdf::pdf(const Intersection &its, const Vector3f &wo, const Vector3f &wi, BxdfType flags) const
{
    if (bxdfs_count == 0)
        return 0;

    Vector3f wo_local = world_to_local(its, wo);
    Vector3f wi_local = world_to_local(its, wi);

    float pdf = 0;
    int count = 0;
    for (int i = 0; i < bxdfs_count; ++i)
    {
        if (bxdfs[i]->matches(flags))
        {
            pdf += bxdfs[i]->pdf(wo_local, wi_local);
            ++count;
        }
    }

    float v = (count > 0) ? pdf / count : 0;
    return v;
}
#endif

Bsdf *make_matte_material(const Spectrum &Kd)
{
    Bsdf *bsdf = new Bsdf();

    if (!Kd.is_black())
        bsdf->add(new LambertianReflection(Kd));

    return bsdf;
}

Bsdf *make_plastic_material(const Spectrum &Kd, const Spectrum &Ks, float roughness)
{
    Bsdf *bsdf = new Bsdf();

    if (!Kd.is_black())
        bsdf->add(new LambertianReflection(Kd));

    if (!Ks.is_black())
    {
        MicrofacetDistribution *distribution = new TrowbridgeReitzDistribution(roughness, roughness);
        Fresnel *fresnel = new FresnelDielectric(1, 1.5);
        MicrofacetReflection *specular = new MicrofacetReflection(Ks, distribution, fresnel);
        bsdf->add(specular);
    }

    return bsdf;
}

Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const
{
    return R * INV_PI;
}

float MicrofacetDistribution::G(const Vector3f &wo, const Vector3f &wi) const
{
    return (float)1 / (1 + lambda(wo) + lambda(wi));
}

float TrowbridgeReitzDistribution::D(const Vector3f &wh) const
{
    float tan2_theta = ::tan2_theta(wh);
    if (std::isinf(tan2_theta))
        return 0;

    float cos4_theta = cos2_theta(wh) * cos2_theta(wh);
    float e = (cos2_phi(wh) / (alphax * alphax) + sin2_phi(wh) / (alphay * alphay)) * tan2_theta;

    return (float)1 / (PI * alphax * alphay * cos4_theta * (1 + e) * (1 + e));
}

float TrowbridgeReitzDistribution::lambda(const Vector3f &w) const
{
    float abs_tan_theta = std::abs(tan_theta(w));
    if (std::isinf(abs_tan_theta))
        return 0;

    float alpha = std::sqrt(cos2_phi(w) * alphax * alphax +
                            sin2_phi(w) * alphay * alphay);

    float alpha2_tan2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
    return (-1 + std::sqrt((float)1 / alpha2_tan2_theta)) / 2;
}

static float fr_dielectric(float cos_theta_i, float eta_i, float eta_t)
{
    cos_theta_i = clamp(cos_theta_i, -1, 1);

    if (cos_theta_i > 0)
    {
        std::swap(eta_i, eta_t);
        cos_theta_i = std::abs(cos_theta_i);
    }

    float sin_theta_i = std::sqrt(std::max((float)0, 1 - cos_theta_i * cos_theta_i));
    float sin_theta_t = eta_i / eta_t * sin_theta_i;
    if (sin_theta_t >= 1)
        return 1;
    float cos_theta_t = std::sqrt(std::max((float)0, 1 - sin_theta_t * sin_theta_t));

    float R_parl = ((eta_t * cos_theta_i) - (eta_i * cos_theta_t)) /
                  ((eta_t * cos_theta_i) + (eta_i * cos_theta_t));
    float R_perp = ((eta_i * cos_theta_i) - (eta_t * cos_theta_t)) /
                  ((eta_i * cos_theta_i) + (eta_t * cos_theta_t));

    return (R_parl * R_parl + R_perp * R_perp) / 2;
}

Spectrum FresnelDielectric::evaluate(float cos_theta_i) const
{
    return fr_dielectric(cos_theta_i, eta_i, eta_t);
}

Spectrum MicrofacetReflection::f(const Vector3f &wo, const Vector3f &wi) const
{
    float cos_theta_0 = abs_cos_theta(wo);
    float cos_theta_i = abs_cos_theta(wi);

    if ((cos_theta_i == 0) || (cos_theta_0 == 0))
        return Spectrum(0);

    Vector3f wh = wi + wo;

    if ((wh.x == 0) && (wh.y == 0) && (wh.z == 0))
        return Spectrum(0);

    wh = normalize(wh);

    Spectrum F = fresnel->evaluate(dot(wi, wh));

    return R * distribution->D(wh) * distribution->G(wo, wi) * F / (4 * cos_theta_i * cos_theta_0);
}
