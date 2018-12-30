#pragma once

#include "spectrum.h"
#include "bsdf.h"

struct Scene;
struct Ray;
struct Camera;
struct Intersection;

struct Integrator
{
    virtual Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0) const = 0;
    virtual void render(const Scene &scene, const Camera &camera) const = 0;

    Spectrum specular_common(const Scene &scene, const Intersection &its, const Ray &ray, int depth, BxdfType type) const;
    Spectrum specular_reflect(const Scene &scene, const Intersection &its, const Ray &ray, int depth) const;
    Spectrum specular_transmit(const Scene &scene, const Intersection &its, const Ray &ray, int depth) const;
};

struct WhittedIntegrator : public Integrator
{
    int max_depth;

    WhittedIntegrator(int max_depth) : max_depth(max_depth) {}

    Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0) const override;
    void render(const Scene &scene, const Camera &camera) const override;
};
