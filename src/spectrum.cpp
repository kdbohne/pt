#include "spectrum.h"

RgbSpectrum::RgbSpectrum(float v)
    : CoefficientSpectrum<3>(v)
{
}

RgbSpectrum::RgbSpectrum(const CoefficientSpectrum<3> &v, SpectrumType type)
    : CoefficientSpectrum<3>(v)
{
}

RgbSpectrum::RgbSpectrum(float r, float g, float b)
{
    c[0] = r;
    c[1] = g;
    c[2] = b;
}

void RgbSpectrum::to_xyz(float xyz[3]) const
{
    rgb_to_xyz(c, xyz);
}

float RgbSpectrum::y() const
{
    const float y_weight[3] = {0.212671, 0.715160, 0.072169};
    return y_weight[0] * c[0] + y_weight[1] * c[1] + y_weight[2] * c[2];
}

void xyz_to_rgb(const float xyz[3], float rgb[3])
{
    rgb[0] =  3.240479 * xyz[0] - 1.537150 * xyz[1] - 0.498535 * xyz[2];
    rgb[1] = -0.969256 * xyz[0] + 1.875991 * xyz[1] + 0.041556 * xyz[2];
    rgb[2] =  0.055648 * xyz[0] - 0.204043 * xyz[1] + 1.057311 * xyz[2];
}

void rgb_to_xyz(const float rgb[3], float xyz[3])
{
    xyz[0] = 0.412453 * rgb[0] + 0.357580 * rgb[1] + 0.180423 * rgb[2];
    xyz[1] = 0.212671 * rgb[0] + 0.715160 * rgb[1] + 0.072169 * rgb[2];
    xyz[2] = 0.019334 * rgb[0] + 0.119193 * rgb[1] + 0.950227 * rgb[2];
}
