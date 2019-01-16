#pragma once

#include "spectrum.h"
#include "bsdf.h"

struct Sampler;
struct Scene;
struct Ray;
class MemoryArena;
struct Camera;
struct Intersection;

struct Integrator
{
    Sampler *sampler;

    virtual Spectrum Li(const Scene &scene, const Ray &ray, MemoryArena &arena, int depth = 0) const = 0;
    virtual void render(const Scene &scene, const Camera &camera) const = 0;

    Spectrum specular_common(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth, BxdfType type) const;
    Spectrum specular_reflect(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth) const;
    Spectrum specular_transmit(const Scene &scene, const Intersection &its, const Ray &ray, MemoryArena &arena, int depth) const;
};

struct WhittedIntegrator : public Integrator
{
    int max_depth;
    Sampler *sampler;

    WhittedIntegrator(int max_depth, Sampler *sampler) : max_depth(max_depth), sampler(sampler) {}

    Spectrum Li(const Scene &scene, const Ray &ray, MemoryArena &arena, int depth = 0) const override;
    void render(const Scene &scene, const Camera &camera) const override;
};
