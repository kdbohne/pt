#pragma once

#include "vector.h"
#include "intersection.h"

#include <memory>

struct Geometry;

struct Light
{
    float Lemit;
    std::shared_ptr<Geometry> geometry;

    Light(float Lemit, const std::shared_ptr<Geometry> &geometry) : Lemit(Lemit), geometry(geometry) {}

    // TODO: return spectrum
    Vector3f L(const Intersection &intersection, const Vector3f &w) const;

#if 0
    // TODO: return spectrum
    // TODO: visibility tester
    Vector3f sample_Li(const Intersection &ref, const Vector2f &u, Vector3f *wi, float *pdf) const;
#endif
};
