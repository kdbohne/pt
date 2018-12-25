#pragma once

#include "spectrum.h"

struct Scene;
struct Ray;
struct Intersection;

struct Integrator
{
    virtual Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0) const = 0;

    Spectrum specular_reflect(const Scene &scene, const Intersection &its, const Ray &ray, int depth) const;
};

struct WhittedIntegrator : public Integrator
{
    int max_depth;

    WhittedIntegrator(int max_depth) : max_depth(max_depth) {}

    Spectrum Li(const Scene &scene, const Ray &ray, int depth = 0) const override;
};
