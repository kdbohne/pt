#include "film.h"
#include <fstream>

#include <OpenEXR/ImfRgba.h>
#include <OpenEXR/ImfRgbaFile.h>

Film::Film(const Point2i &resolution)
    : resolution(resolution), pixels(resolution.x * resolution.y)
{
}

void Film::set_pixel(int x, int y, const Spectrum &v)
{
    assert((x >= 0) && (x < resolution.x));
    assert((y >= 0) && (y < resolution.y));

    float xyz[3];
    v.to_xyz(xyz);

    for (int i = 0; i < 3; ++i)
        pixels[y * resolution.x + x].xyz[i] = xyz[i];
}

static void tone_map_reinhard(float *rgb)
{
    rgb[0] /= 1 + rgb[0];
    rgb[1] /= 1 + rgb[1];
    rgb[2] /= 1 + rgb[2];
}

static void accurate_linear_to_srgb(float *linear)
{
    linear[0] = (linear[0] <= 0.0031308) ? linear[0] * 12.92 : (std::pow(std::abs(linear[0]), 1.0 / 2.4) * 1.055) - 0.055;
    linear[1] = (linear[1] <= 0.0031308) ? linear[1] * 12.92 : (std::pow(std::abs(linear[1]), 1.0 / 2.4) * 1.055) - 0.055;
    linear[2] = (linear[2] <= 0.0031308) ? linear[2] * 12.92 : (std::pow(std::abs(linear[2]), 1.0 / 2.4) * 1.055) - 0.055;
}

void Film::write_ppm(const std::string &path) const
{
    std::ofstream file(path);

    file << "P3\n";
    file << resolution.x << " " << resolution.y << "\n";
    file << "255" << "\n";

    for (int y = 0; y < resolution.y; ++y)
    {
        for (int x = 0; x < resolution.x; ++x)
        {
            const Pixel *p = &pixels[y * resolution.x + x];

            float rgb[3];
            xyz_to_rgb(p->xyz, rgb);

            tone_map_reinhard(rgb);

            accurate_linear_to_srgb(rgb);

            int r = (int)(rgb[0] * 255);
            int g = (int)(rgb[1] * 255);
            int b = (int)(rgb[2] * 255);

            file << r << " " << g << " " << b << " ";
        }

        file << "\n";
    }
}

void Film::write_exr(const std::string &path) const
{
    Imf::Rgba *rgba = new Imf::Rgba[resolution.x * resolution.y];
    for (int i = 0; i < resolution.x * resolution.y; ++i)
    {
        float rgb[3];
        xyz_to_rgb(pixels[i].xyz, rgb);

        rgba[i] = Imf::Rgba(rgb[0], rgb[1], rgb[2]);
    }

    Imf::RgbaOutputFile file(path.c_str(), resolution.x, resolution.y, Imf::WRITE_RGB);
    file.setFrameBuffer(rgba, 1, resolution.x);
    file.writePixels(resolution.y);
}
