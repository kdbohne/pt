#pragma once

#include "point.h"
#include "spectrum.h"

#include <vector>
#include <string>

struct Pixel
{
    float xyz[3];
};

struct Film
{
    Point2i resolution;
    std::vector<Pixel> pixels;

    Film(const Point2i &resolution);

    void set_pixel(int x, int y, const Spectrum &v);

    void write_ppm(const std::string &path) const;
};
