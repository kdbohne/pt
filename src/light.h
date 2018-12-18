#pragma once

#include "vector.h"
#include "intersection.h"
#include "spectrum.h"

#include <memory>

struct Geometry;

struct Light
{
    Spectrum Lemit;
    std::shared_ptr<Geometry> geometry;

    Light(const Spectrum &Lemit, const std::shared_ptr<Geometry> &geometry);

    Spectrum L(const Intersection &intersection, const Vector3f &w) const;

#if 0
    // TODO: visibility tester
    Spectrum sample_Li(const Intersection &ref, const Vector2f &u, Vector3f *wi, float *pdf) const;
#endif
};
