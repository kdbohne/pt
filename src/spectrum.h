#pragma once

#include "common.h"
#include "math.h"

template<int N>
struct CoefficientSpectrum
{
    float c[N];

    CoefficientSpectrum(float v = 0.0)
    {
        for (int i = 0; i < N; ++i)
            c[i] = v;
    }

    CoefficientSpectrum<N> operator+(const CoefficientSpectrum<N> &s) const
    {
        CoefficientSpectrum<N> result = *this;
        for (int i = 0; i < N; ++i)
            result.c[i] += s.c[i];

        return result;
    }

    CoefficientSpectrum<N> &operator+=(const CoefficientSpectrum<N> &s)
    {
        for (int i = 0; i < N; ++i)
            c[i] += s.c[i];

        return *this;
    }

    CoefficientSpectrum<N> operator-(const CoefficientSpectrum<N> &s) const
    {
        CoefficientSpectrum<N> result = *this;
        for (int i = 0; i < N; ++i)
            result.c[i] -= s.c[i];

        return result;
    }

    CoefficientSpectrum<N> operator*(const CoefficientSpectrum<N> &s) const
    {
        CoefficientSpectrum<N> result = *this;
        for (int i = 0; i < N; ++i)
            result.c[i] *= s.c[i];

        return result;
    }

    CoefficientSpectrum<N> operator*(float s) const
    {
        CoefficientSpectrum<N> result = *this;
        for (int i = 0; i < N; ++i)
            result.c[i] *= s;

        return result;
    }

    CoefficientSpectrum<N> &operator*=(float s)
    {
        for (int i = 0; i < N; ++i)
            c[i] *= s;

        return *this;
    }

    CoefficientSpectrum<N> operator/(float s) const
    {
        assert(s != 0);

        CoefficientSpectrum<N> result = *this;
        for (int i = 0; i < N; ++i)
            result.c[i] /= s;

        return result;
    }

    CoefficientSpectrum<N> &operator/=(float s)
    {
        assert(s != 0);

        for (int i = 0; i < N; ++i)
            c[i] /= s;

        return *this;
    }

    bool is_black() const
    {
        for (int i = 0; i < N; ++i)
        {
            if (c[i] != 0)
                return false;
        }
        return true;
    }

    CoefficientSpectrum<N> clamp(float min = 0, float max = INFINITY) const
    {
        CoefficientSpectrum<N> result;
        for (int i = 0; i < N; ++i)
            result.c[i] = ::clamp(c[i], min, max);

        return result;
    }
};

template<int N>
inline CoefficientSpectrum<N> operator*(float s, const CoefficientSpectrum<N> &cs)
{
    return cs * s;
}

struct RgbSpectrum : public CoefficientSpectrum<3>
{
    RgbSpectrum(float v = 0.0);
    RgbSpectrum(float r, float g, float b);
    RgbSpectrum(const CoefficientSpectrum<3> &v);

    void to_xyz(float xyz[3]) const;
};

// TODO: SampledSpectrum

typedef RgbSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;

void xyz_to_rgb(const float xyz[3], float rgb[3]);
void rgb_to_xyz(const float rgb[3], float xyz[3]);
