#include "bsdf.h"
#include "math.h"
#include "sampling.h"
#include "intersection.h"

bool Bxdf::matches(BxdfType t) const
{
    return (type & t) == type;
}

Spectrum Bxdf::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
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
        if (sampled_type)
            *sampled_type = (BxdfType)0;
        return Spectrum(0);
    }

    int component = std::min((int)std::floor(u[0] * components_count), components_count - 1);

    // TODO CLEANUP: ugly iteration
    Bxdf *bxdf = nullptr;
    int count = component;
    for (int i = 0; i < bxdfs_count; ++i)
    {
        if (bxdfs[i]->matches(type) && (count-- == 0))
        {
            bxdf = bxdfs[i];
            break;
        }
    }
    assert(bxdf);

    Point2f u_remapped(std::min(u[0] * components_count - component, ONE_MINUS_EPSILON), u[1]);

    if (sampled_type)
        *sampled_type = bxdf->type;

    Vector3f wo_local = world_to_local(its, wo);
    if (wo_local.z == 0)
        return 0;

    *pdf = 0;
    Vector3f wi_local;
    Spectrum f = bxdf->sample_f(its, wo_local, &wi_local, u_remapped, pdf, sampled_type);
    if (*pdf == 0)
    {
        if (sampled_type)
            *sampled_type = (BxdfType)0;
        return 0;
    }

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

    if (!(bxdf->type & BSDF_SPECULAR))
    {
        bool reflect = (dot(*wi, its.n) * dot(wo, its.n) > 0);

        f = 0;
        for (int i = 0; i < bxdfs_count; ++i)
        {
            if (bxdfs[i]->matches(type))
            {
                if ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                    (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION)))
                {
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

Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const
{
    return R * INV_PI;
}

float MicrofacetDistribution::G1(const Vector3f &w) const
{
    return 1 / (1 + lambda(w));
}

float MicrofacetDistribution::G(const Vector3f &wo, const Vector3f &wi) const
{
    return (float)1 / (1 + lambda(wo) + lambda(wi));
}

float TrowbridgeReitzDistribution::roughness_to_alpha(float roughness)
{
    roughness = std::max(roughness, (float)1e-3);
    float x = std::log(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x + 0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}

float TrowbridgeReitzDistribution::D(const Vector3f &wh) const
{
    float tan2_theta = ::tan2_theta(wh);
    if (std::isinf(tan2_theta))
        return 0;

    float cos4_theta = cos2_theta(wh) * cos2_theta(wh);
    float e = (cos2_phi(wh) / (alpha_x * alpha_x) + sin2_phi(wh) / (alpha_y * alpha_y)) * tan2_theta;

    return (float)1 / (PI * alpha_x * alpha_y * cos4_theta * (1 + e) * (1 + e));
}

float TrowbridgeReitzDistribution::lambda(const Vector3f &w) const
{
    float abs_tan_theta = std::abs(tan_theta(w));
    if (std::isinf(abs_tan_theta))
        return 0;

    float alpha = std::sqrt(cos2_phi(w) * alpha_x * alpha_x +
                            sin2_phi(w) * alpha_y * alpha_y);

    float alpha2_tan2_theta = (alpha * abs_tan_theta) * (alpha * abs_tan_theta);
    return (-1 + std::sqrt((float)1 / alpha2_tan2_theta)) / 2;
}

static void trowbridge_reitz_sample_11(float cos_theta, float u0, float u1, float *slope_x, float *slope_y)
{
    if (cos_theta > 0.9999)
    {
        float r = std::sqrt(u0 / (1 - u0));
        float phi = 6.28318530718 * u1;

        *slope_x = r * cos(phi);
        *slope_y = r * sin(phi);

        return;
    }

    float sin_theta = std::sqrt(std::max((float)0, (float)1 - cos_theta * cos_theta));
    float tan_theta = sin_theta / cos_theta;
    float a = 1 / tan_theta;
    float G1 = 2 / (1 + std::sqrt(1 + 1 / (a * a)));

    float A = 2 * u0 / G1 - 1;
    float tmp = 1 / (A * A - 1);
    if (tmp > 1e10)
        tmp = 1e10;
    float B = tan_theta;
    float D = std::sqrt(std::max((float)(B * B * tmp * tmp - (A * A - B * B) * tmp), (float)0));

    float slope_x_1 = B * tmp - D;
    float slope_x_2 = B * tmp + D;

    *slope_x = ((A < 0) || (slope_x_2 > 1 / tan_theta)) ? slope_x_1 : slope_x_2;

    float S;
    if (u1 > 0.5f)
    {
        S = 1;
        u1 = 2 * (u1 - 0.5);
    }
    else
    {
        S = -1;
        u1 = 2 * (0.5 - u1);
    }

    float z = (u1 * (u1 * (u1 * 0.27385 - 0.73369) + 0.46341)) /
              (u1 * (u1 * (u1 * 0.093073 + 0.309420) - 1.000000) + 0.597999);

    *slope_y = S * z * std::sqrt(1 + *slope_x * *slope_x);

    assert(!std::isinf(*slope_y));
    assert(!std::isnan(*slope_y));
}

static Vector3f trowbridge_reitz_sample(const Vector3f &wi, float alpha_x, float alpha_y, float u0, float u1)
{
    Vector3f wi_stretched = normalize(Vector3f(alpha_x * wi.x, alpha_y * wi.y, wi.z));

    float slope_x, slope_y;
    trowbridge_reitz_sample_11(cos_theta(wi_stretched), u0, u1, &slope_x, &slope_y);

    float tmp = cos_phi(wi_stretched) * slope_x - sin_phi(wi_stretched) * slope_y;
    slope_y = sin_phi(wi_stretched) * slope_x + cos_phi(wi_stretched) * slope_y;
    slope_x = tmp;

    slope_x = alpha_x * slope_x;
    slope_y = alpha_y * slope_y;

    return normalize(Vector3f(-slope_x, -slope_y, 1));
}

Vector3f TrowbridgeReitzDistribution::sample_wh(const Vector3f &wo, const Point2f &u) const
{
    Vector3f wh;
    if (!sample_visible_area)
    {
        float cos_theta = 0;
        float phi = (2 * PI) * u[1];
        if (alpha_x == alpha_y)
        {
            float tan_theta2 = alpha_x * alpha_x * u[0] / (1 - u[0]);
            cos_theta = 1 / std::sqrt(1 + tan_theta2);
        }
        else
        {
            phi = std::atan(alpha_y / alpha_x * std::tan(2 * PI * u[1] + 0.5 * PI));
            if (u[1] > 0.5)
                phi += PI;

            float sin_phi = std::sin(phi);
            float cos_phi = std::cos(phi);

            float alpha_x2 = alpha_x * alpha_x;
            float alpha_y2 = alpha_y * alpha_y;
            float alpha2 = 1 / (cos_phi * cos_phi / alpha_x2 + sin_phi * sin_phi / alpha_y2);

            float tan_theta2 = alpha2 * u[0] / (1 - u[0]);
            cos_theta = 1 / std::sqrt(1 + tan_theta2);
        }

        float sin_theta = std::sqrt(std::max((float)0, (float)1 - cos_theta * cos_theta));

        wh = spherical_direction(sin_theta, cos_theta, phi);
        if (!same_hemisphere(wo, wh))
            wh = -wh;
    }
    else
    {
        bool flip = (wo.z < 0);

        wh = trowbridge_reitz_sample(flip ? -wo : wo, alpha_x, alpha_y, u[0], u[1]);
        if (flip)
            wh = -wh;
    }

    return wh;
}

float MicrofacetDistribution::pdf(const Vector3f &wo, const Vector3f &wh) const
{
    if (sample_visible_area)
        return D(wh) * G1(wo) * abs_dot(wo, wh) / abs_cos_theta(wo);
    else
        return D(wh) * abs_cos_theta(wh);
}

static float fr_dielectric(float cos_theta_i, float eta_i, float eta_t)
{
    cos_theta_i = clamp(cos_theta_i, -1, 1);

    bool entering = (cos_theta_i > 0);
    if (!entering)
    {
        std::swap(eta_i, eta_t);
        cos_theta_i = std::abs(cos_theta_i);
    }

    float sin_theta_i = std::sqrt(std::max((float)0, (float)1 - cos_theta_i * cos_theta_i));
    float sin_theta_t = eta_i / eta_t * sin_theta_i;
    if (sin_theta_t >= 1)
        return 1;
    float cos_theta_t = std::sqrt(std::max((float)0, (float)1 - sin_theta_t * sin_theta_t));

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

Spectrum FresnelSpecular::f(const Vector3f &wo, const Vector3f &wi) const
{
    return Spectrum(0);
}

Spectrum FresnelSpecular::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
{
    float F = fr_dielectric(cos_theta(wo), eta_a, eta_b);
    if (u[0] < F)
    {
        // Specular reflection.

        // Perfect specular reflection direction.
        *wi = Vector3f(-wo.x, -wo.y, wo.z);
        *pdf = F;

        if (sampled_type)
            *sampled_type = (BxdfType)(BSDF_SPECULAR | BSDF_REFLECTION);

        return F * R / abs_cos_theta(*wi);
    }
    else
    {
        // Specular transmission.

        // Determine incident and transmitted etas.
        bool entering = (cos_theta(wo) > 0);
        float eta_i = entering ? eta_a : eta_b;
        float eta_t = entering ? eta_b : eta_a;

        if (!refract(wo, face_forward(Normal3f(0, 0, 1), wo), eta_i / eta_t, wi))
            return 0;

        Spectrum ft = T * (1 - F);

        // TODO: actually check transport mode here instead of this always-true condition
        if (true)
            ft *= (eta_i * eta_i) / (eta_t * eta_t);

        if (sampled_type)
            *sampled_type = (BxdfType)(BSDF_SPECULAR | BSDF_TRANSMISSION);

        *pdf = 1 - F;

        return ft / abs_cos_theta(*wi);
    }
}

float FresnelSpecular::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    return 0;
}

Spectrum SpecularReflection::f(const Vector3f &wo, const Vector3f &wi) const
{
    return Spectrum(0);
}

Spectrum SpecularReflection::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
{
    // Perfect specular reflection direction.
    *wi = Vector3f(-wo.x, -wo.y, wo.z);
    *pdf = 1;

    return fresnel->evaluate(cos_theta(*wi)) * R / abs_cos_theta(*wi);
}

float SpecularReflection::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    return 0;
}

Spectrum SpecularTransmission::f(const Vector3f &wo, const Vector3f &wi) const
{
    return Spectrum(0);
}

Spectrum SpecularTransmission::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
{
    // Determine incident and transmitted etas.
    bool entering = (cos_theta(wo) > 0);
    float eta_i = entering ? eta_a : eta_b;
    float eta_t = entering ? eta_b : eta_a;

    if (!refract(wo, face_forward(Normal3f(0, 0, 1), wo), eta_i / eta_t, wi))
        return 0;

    *pdf = 1;

    Spectrum ft = T * (Spectrum(1) - fresnel.evaluate(cos_theta(*wi)));

    // TODO: actually check transport mode here instead of this always-true condition
    if (true)
        ft *= (eta_i * eta_i) / (eta_t * eta_t);

    return ft / abs_cos_theta(*wi);
}

float SpecularTransmission::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    return 0;
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

Spectrum MicrofacetReflection::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
{
    if (wo.z == 0)
        return 0;

    Vector3f wh = distribution->sample_wh(wo, u);
    *wi = reflect(wo, wh);
    if (!same_hemisphere(wo, *wi))
        return Spectrum(0);

    *pdf = distribution->pdf(wo, wh) / (4 * dot(wo, wh));

    return f(wo, *wi);
}

float MicrofacetReflection::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    if (!same_hemisphere(wo, wi))
        return 0;

    Vector3f wh = normalize(wo + wi);
    return distribution->pdf(wo, wh) / (4 * dot(wo, wh));
}

Spectrum MicrofacetTransmission::f(const Vector3f &wo, const Vector3f &wi) const
{
    if (same_hemisphere(wo, wi))
        return Spectrum(0);

    float cos_theta_o = cos_theta(wo);
    float cos_theta_i = cos_theta(wi);
    if ((cos_theta_o == 0) || (cos_theta_i == 0))
        return Spectrum(0);

    float eta = (cos_theta_o > 0) ? (eta_b / eta_a) : (eta_a / eta_b);
    Vector3f wh = normalize(wo + wi * eta);
    if (wh.z < 0)
        wh = -wh;

    Spectrum F = fresnel.evaluate(dot(wo, wh));

    float sqrt_denom = dot(wo, wh) + eta * dot(wi, wh);
    // TODO: actually check transport mode here instead of this always-true condition
    float factor = (true) ? (1 / eta) : 1;

    return (Spectrum(1) - F) * T *
           std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
           abs_dot(wi, wh) * abs_dot(wo, wh) * factor * factor /
           (cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom));
}

Spectrum MicrofacetTransmission::sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType *sampled_type) const
{
    if (wo.z == 0)
        return 0;

    Vector3f wh = distribution->sample_wh(wo, u);
    float eta = (cos_theta(wo) > 0) ? (eta_a / eta_b) : (eta_b / eta_a);
    if (!refract(wo, Normal3f(wh), eta, wi))
        return 0;

    *pdf = this->pdf(wo, *wi);

    return f(wo, *wi);
}

float MicrofacetTransmission::pdf(const Vector3f &wo, const Vector3f &wi) const
{
    if (!same_hemisphere(wo, wi))
        return 0;

    float eta = (cos_theta(wo) > 0) ? (eta_b / eta_a) : (eta_a / eta_b);
    Vector3f wh = normalize(wo + wi * eta);

    float sqrt_denom = dot(wo, wh) + eta * dot(wi, wh);
    float dwh_dwi = std::abs((eta * eta * dot(wi, wh)) / (sqrt_denom * sqrt_denom));

    return distribution->pdf(wo, wh) / dwh_dwi;
}
