#pragma once

#include "vector.h"
#include <vector>
#include <string>

struct Pixel
{
    // TODO: spectrum
    Vector3f rgb;
};

struct Film
{
    Vector2i resolution;
    std::vector<Pixel> pixels;

    Film(Vector2i resolution) : resolution(resolution), pixels(resolution.x * resolution.y) {}

    void write_ppm(const std::string &path) const;
};
