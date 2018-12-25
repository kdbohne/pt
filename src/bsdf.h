#pragma once

#include "spectrum.h"
#include "vector.h"
#include "point.h"
#include "math.h"

struct Intersection;

enum BxdfType
{
    // NOTE: only handling reflection for now
    // TODO: transmission

    // Must be one of the following base types:
    BSDF_REFLECTION   = 0x1,
    BSDF_TRANSMISSION = 0x2,

    // Must be one of the following types:
    BSDF_DIFFUSE  = 0x4,
    BSDF_GLOSSY   = 0x8,
    BSDF_SPECULAR = 0x10,

    BSDF_ALL = BSDF_REFLECTION | BSDF_TRANSMISSION |
               BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
};

struct Bxdf
{
    BxdfType type;

    Bxdf(BxdfType type) : type(type) {}

    bool matches(BxdfType t) const;

    virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;
    virtual Spectrum sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf) const;
    virtual float pdf(const Vector3f &wo, const Vector3f &wi) const;
};

struct Bsdf
{
    static constexpr int MAX_BXDFS_COUNT = 8; // TODO: ?
    Bxdf *bxdfs[MAX_BXDFS_COUNT];
    int bxdfs_count;

    Bsdf() : bxdfs_count(0) {}

    void add(Bxdf *bxdf);

    int get_matching_components_count(BxdfType flags) const;

    Spectrum f(const Intersection &its, const Vector3f &wo, const Vector3f &wi, BxdfType flags = BSDF_ALL) const;
    Spectrum sample_f(const Intersection &its, const Vector3f &wo, Vector3f *wi, const Point2f &u, float *pdf, BxdfType type = BSDF_ALL, BxdfType *sampled_type = nullptr) const;
};

struct LambertianReflection : public Bxdf
{
    const Spectrum R;

    LambertianReflection(const Spectrum &R) : Bxdf((BxdfType)(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const override;
};

struct MicrofacetDistribution
{
    bool sample_visible_area;

    MicrofacetDistribution(bool sample_visible_area) : sample_visible_area(sample_visible_area) {}

    virtual float D(const Vector3f &wh) const = 0;
    virtual float lambda(const Vector3f &w) const = 0;

    float G(const Vector3f &wo, const Vector3f &wi) const;
};

struct TrowbridgeReitzDistribution : public MicrofacetDistribution
{
    float alphax, alphay;

    TrowbridgeReitzDistribution(float ax, float ay, bool sample_visible_area = true)
        : MicrofacetDistribution(sample_visible_area), alphax(ax), alphay(ay) {}

    float D(const Vector3f &wh) const override;
    float lambda(const Vector3f &w) const override;
};

struct Fresnel
{
    virtual Spectrum evaluate(float cos_i) const = 0;
};

struct FresnelDielectric : public Fresnel
{
    float eta_i, eta_t;

    FresnelDielectric(float eta_i, float eta_t) : eta_i(eta_i), eta_t(eta_t) {}

    Spectrum evaluate(float cos_theta_i) const override;
};

struct MicrofacetReflection : public Bxdf
{
    const Spectrum R;
    const MicrofacetDistribution *distribution;
    const Fresnel *fresnel;

    MicrofacetReflection(const Spectrum &R, const MicrofacetDistribution *distribution, const Fresnel *fresnel)
        : Bxdf((BxdfType)(BSDF_REFLECTION | BSDF_GLOSSY)),
          R(R), distribution(distribution), fresnel(fresnel) {}

    Spectrum f(const Vector3f &wo, const Vector3f &wi) const override;
};

inline float cos_theta(const Vector3f &w)
{
    return w.z;
}

inline float cos2_theta(const Vector3f &w)
{
    return w.z * w.z;
}

inline float abs_cos_theta(const Vector3f &w)
{
    return std::abs(w.z);
}

inline float sin2_theta(const Vector3f &w)
{
    return std::max((float)0, (float)1 - cos2_theta(w));
}

inline float sin_theta(const Vector3f &w)
{
    return std::sqrt(sin2_theta(w));
}

inline float tan_theta(const Vector3f &w)
{
    return sin_theta(w) / cos_theta(w);
}

inline float tan2_theta(const Vector3f &w)
{
    return sin2_theta(w) / cos2_theta(w);
}

inline float cos_phi(const Vector3f &w)
{
    float sin_theta = ::sin_theta(w);
    if (sin_theta == 0)
        return 1;

    return clamp(w.x / sin_theta, -1, 1);
}

inline float cos2_phi(const Vector3f &w)
{
    return cos_phi(w) * cos_phi(w);
}

inline float sin_phi(const Vector3f &w)
{
    float sin_theta = ::sin_theta(w);
    if (sin_theta == 0)
        return 0;

    return clamp(w.y / sin_theta, -1, 1);
}

inline float sin2_phi(const Vector3f &w)
{
    return sin_phi(w) * sin_phi(w);
}

inline bool same_hemisphere(const Vector3f &w, const Vector3f &wp)
{
    return (w.z * wp.z > 0);
}
