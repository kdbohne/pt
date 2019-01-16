#include "material.h"
#include "bsdf.h"
#include "texture.h"
#include "intersection.h"
#include "memory.h"

// TODO MEMORY: use a pool allocator or memory arena for each allocation here

void MatteMaterial::evaluate_surface(const Ray &ray, MemoryArena &arena, Intersection *its) const
{
    its->bsdf = ARENA_ALLOC(arena, Bsdf)();

    Spectrum r = Kd->evaluate(*its).clamp();
    if (!r.is_black())
        its->bsdf->add(ARENA_ALLOC(arena, LambertianReflection)(r));
}

void PlasticMaterial::evaluate_surface(const Ray &ray, MemoryArena &arena, Intersection *its) const
{
    its->bsdf = ARENA_ALLOC(arena, Bsdf)();

    Spectrum kd = Kd->evaluate(*its).clamp();
    if (!kd.is_black())
        its->bsdf->add(ARENA_ALLOC(arena, LambertianReflection)(kd));

    Spectrum ks = Ks->evaluate(*its).clamp();
    if (!ks.is_black())
    {
        float rough = roughness->evaluate(*its);
        if (remap_roughness)
            rough = TrowbridgeReitzDistribution::roughness_to_alpha(rough);

        MicrofacetDistribution *distribution = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(rough, rough);
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1, 1.5);
        MicrofacetReflection *specular = ARENA_ALLOC(arena, MicrofacetReflection)(ks, distribution, fresnel);
        its->bsdf->add(specular);
    }
}

void GlassMaterial::evaluate_surface(const Ray &ray, MemoryArena &arena, Intersection *its) const
{
    its->bsdf = ARENA_ALLOC(arena, Bsdf)();

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
        its->bsdf->add(ARENA_ALLOC(arena, FresnelSpecular)(kr, kt, 1, eta));
        return;
    }

    if (remap_roughness)
    {
        u_rough = TrowbridgeReitzDistribution::roughness_to_alpha(u_rough);
        v_rough = TrowbridgeReitzDistribution::roughness_to_alpha(v_rough);
    }

    MicrofacetDistribution *distribution = nullptr;
    if (!is_specular)
        distribution = ARENA_ALLOC(arena, TrowbridgeReitzDistribution)(u_rough, v_rough);

    if (!kr.is_black())
    {
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1, eta);
        if (is_specular)
            its->bsdf->add(ARENA_ALLOC(arena, SpecularReflection)(kr, fresnel));
        else
            its->bsdf->add(ARENA_ALLOC(arena, MicrofacetReflection)(kr, distribution, fresnel));
    }

    if (!kt.is_black())
    {
        if (is_specular)
            its->bsdf->add(ARENA_ALLOC(arena, SpecularTransmission)(kt, 1, eta));
        else
            its->bsdf->add(ARENA_ALLOC(arena, MicrofacetTransmission)(kt, distribution, 1, eta));
    }
}
