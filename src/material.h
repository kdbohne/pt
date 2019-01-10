#pragma once

#include "spectrum.h"

struct Ray;
struct Intersection;
template<typename T> struct Texture;

struct Material
{
    virtual void evaluate_surface(const Ray &ray, Intersection *its) const = 0;
};

struct MatteMaterial : public Material
{
    Texture<Spectrum> *Kd;

    MatteMaterial(Texture<Spectrum> *Kd) : Kd(Kd) {}

    void evaluate_surface(const Ray &ray, Intersection *its) const override;
};

struct PlasticMaterial : public Material
{
    Texture<Spectrum> *Kd;
    Texture<Spectrum> *Ks;
    Texture<float> *roughness;
    bool remap_roughness;

    PlasticMaterial(Texture<Spectrum> *Kd, Texture<Spectrum> *Ks, Texture<float> *roughness, bool remap_roughness)
        : Kd(Kd), Ks(Ks), roughness(roughness), remap_roughness(remap_roughness) {}

    void evaluate_surface(const Ray &ray, Intersection *its) const override;
};

struct GlassMaterial : public Material
{
    Texture<Spectrum> *Kr;
    Texture<Spectrum> *Kt;
    Texture<float> *u_roughness;
    Texture<float> *v_roughness;
    Texture<float> *index;
    bool remap_roughness;

    GlassMaterial(Texture<Spectrum> *Kr, Texture<Spectrum> *Kt, Texture<float> *u_roughness, Texture<float> *v_roughness, Texture<float> *index, bool remap_roughness)
        : Kr(Kr), Kt(Kt), u_roughness(u_roughness), v_roughness(v_roughness), index(index), remap_roughness(remap_roughness) {}

    void evaluate_surface(const Ray &ray, Intersection *its) const override;
};
