#pragma once

#include "spectrum.h"

struct Bxdf
{
    virtual Spectrum f(/* TODO */) const = 0;
};

struct LambertianReflection : public Bxdf
{
    Spectrum R;

    LambertianReflection(const Spectrum &R);

    Spectrum f(/* TODO */) const override;
};

struct Bsdf
{
    static constexpr int MAX_BXDFS_COUNT = 8; // TODO: ?
    Bxdf *bxdfs[MAX_BXDFS_COUNT];
    int bxdfs_count;

    Bsdf() : bxdfs_count(0) {}

    void add(Bxdf *bxdf);

    Spectrum f(/* TODO */) const;
};
