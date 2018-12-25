#include "bsdf.h"
#include "math.h"

LambertianReflection::LambertianReflection(const Spectrum &R)
    : R(R)
{
}

Spectrum LambertianReflection::f(/* TODO */) const
{
    return R * INV_PI;
}

void Bsdf::add(Bxdf *bxdf)
{
    assert(bxdfs_count < MAX_BXDFS_COUNT);
    bxdfs[bxdfs_count++] = bxdf;
}

Spectrum Bsdf::f(/* TODO */) const
{
    Spectrum f(0);
    for (int i = 0; i < bxdfs_count; ++i)
        f += bxdfs[i]->f(/* TODO */);

    return f;
}
