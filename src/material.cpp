#include "material.h"
#include "bsdf.h"
#include "texture.h"
#include "intersection.h"

// TODO MEMORY: use a pool allocator or memory arena for each allocation here

void MatteMaterial::evaluate_surface(const Ray &ray, Intersection *its) const
{
    its->bsdf = new Bsdf();

    Spectrum r = Kd->evaluate(*its).clamp();
    if (!r.is_black())
        its->bsdf->add(new LambertianReflection(r));
}

void PlasticMaterial::evaluate_surface(const Ray &ray, Intersection *its) const
{
    its->bsdf = new Bsdf();

    Spectrum kd = Kd->evaluate(*its).clamp();
    if (!kd.is_black())
        its->bsdf->add(new LambertianReflection(kd));

    Spectrum ks = Ks->evaluate(*its).clamp();
    if (!ks.is_black())
    {
        float rough = roughness->evaluate(*its);
        if (remap_roughness)
            rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);

        MicrofacetDistribution *distribution = new TrowbridgeReitzDistribution(rough, rough);
        Fresnel *fresnel = new FresnelDielectric(1, 1.5);
        MicrofacetReflection *specular = new MicrofacetReflection(ks, distribution, fresnel);
        its->bsdf->add(specular);
    }
}

void GlassMaterial::evaluate_surface(const Ray &ray, Intersection *its) const
{
    its->bsdf = new Bsdf();

    Spectrum kr = Kr->evaluate(*its).clamp();
    Spectrum kt = Kt->evaluate(*its).clamp();
    float u_rough = u_roughness->evaluate(*its);
    float v_rough = v_roughness->evaluate(*its);
    float eta = index->evaluate(*its);

    if (kr.is_black() && kt.is_black())
        return;

    bool is_specular = (u_rough == 0) && (v_rough == 0);
    bool allow_multiple_lobes = false; // TODO: configurable?
    if (is_specular && allow_multiple_lobes)
    {
        its->bsdf->add(new FresnelSpecular(kr, kt, 1, eta));
        return;
    }

    if (remap_roughness)
    {
        u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
        v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
    }

    MicrofacetDistribution *distribution = nullptr;
    if (!is_specular)
        distribution = new TrowbridgeReitzDistribution(u_rough, v_rough);

    if (!kr.is_black())
    {
        Fresnel *fresnel = new FresnelDielectric(1, eta);
        if (is_specular)
            its->bsdf->add(new SpecularReflection(kr, fresnel));
        else
            its->bsdf->add(new MicrofacetReflection(kr, distribution, fresnel));
    }

    if (!kt.is_black())
    {
        if (is_specular)
            its->bsdf->add(new SpecularTransmission(kt, 1, eta));
        else
            its->bsdf->add(new MicrofacetTransmission(kt, distribution, 1, eta));
    }
}
