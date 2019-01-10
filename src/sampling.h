#pragma once

#include "vector.h"
#include "point.h"

#include <vector>

float uniform_float();

Vector3f uniform_sample_sphere(const Point2f &u);
Point2f uniform_sample_triangle(const Point2f &u);
Point2f concentric_sample_disk(const Point2f &u);
Vector3f cosine_sample_hemisphere(const Point2f &u);

float uniform_cone_pdf(float cos_theta_max);

struct Distribution1d
{
    std::vector<float> func, cdf;
    float func_int;

    Distribution1d(const float *f, int n);

    float sample_continuous(float u, float *pdf, int *offset = nullptr) const;
};

struct Distribution2d
{
    std::vector<Distribution1d *> p_conditional_v;
    Distribution1d *p_marginal;

    Distribution2d(float *func, int nu, int nv);

    Point2f sample_continuous(const Point2f &u, float *pdf) const;
};
