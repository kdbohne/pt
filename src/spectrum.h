#pragma once

template<int N>
struct CoefficientSpectrum
{
    static const int samples_count = N;
    float c[N];

    CoefficientSpectrum(float v = 0.0);

    virtual void to_xyz(float xyz[3]) const = 0;
};

struct RgbSpectrum : public CoefficientSpectrum<3>
{
    RgbSpectrum(float v = 0.0);
    RgbSpectrum(const CoefficientSpectrum<3> &v);

    void to_xyz(float *xyz) const override;
};

// TODO: SampledSpectrum

typedef RgbSpectrum Spectrum;
//typedef SampledSpectrum Spectrum;

void xyz_to_rgb(const float xyz[3], float rgb[3]);
void rgb_to_xyz(const float rgb[3], float xyz[3]);
